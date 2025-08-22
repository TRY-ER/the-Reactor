"""
Simple Web Scraper - Master Function for Single URL Scraping (Synchronous Version)
"""

import os
import re
import json
from pathlib import Path
from urllib.parse import urljoin, urlparse, unquote
from typing import Dict, List, Optional
from datetime import datetime
import hashlib

import requests
from bs4 import BeautifulSoup
from PIL import Image
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def sanitize_filename(filename: str) -> str:
    """Sanitize filename for filesystem compatibility"""
    filename = re.sub(r'[^\w\s.-]', '_', filename)
    filename = re.sub(r'\s+', '_', filename)
    filename = filename.strip('._')
    
    if len(filename) > 100:
        filename = filename[:100]
    
    return filename or "unnamed"


def create_unique_output_dir(base_url: str, base_dir: str = "scraped_data") -> Path:
    """Create a unique output directory based on URL and timestamp"""
    # Extract domain from URL
    parsed_url = urlparse(base_url)
    domain = parsed_url.netloc.replace('www.', '')
    
    # Create timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Create unique hash from URL
    url_hash = hashlib.md5(base_url.encode()).hexdigest()[:8]
    
    # Create directory name
    dir_name = f"{sanitize_filename(domain)}_{timestamp}_{url_hash}"
    output_dir = Path(base_dir) / dir_name
    
    # Create directories
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "images").mkdir(exist_ok=True)
    
    return output_dir


def download_image(image_url: str, base_url: str, output_dir: Path, session: requests.Session) -> Optional[str]:
    """Download a single image"""
    try:
        # Resolve relative URLs
        full_url = urljoin(base_url, image_url)
        
        # Get filename
        parsed_url = urlparse(full_url)
        filename = os.path.basename(unquote(parsed_url.path))
        
        if not filename or '.' not in filename:
            filename = f"image_{hash(full_url) % 1000000}.jpg"
        
        # Check file extension
        allowed_formats = {'.jpg', '.jpeg', '.png', '.gif', '.webp', '.svg'}
        file_ext = Path(filename).suffix.lower()
        if file_ext not in allowed_formats:
            return None
        
        file_path = output_dir / "images" / sanitize_filename(filename)
        
        # Skip if already exists
        if file_path.exists():
            return str(file_path)
        
        # Download image
        response = session.get(full_url, timeout=10)
        if response.status_code == 200:
            content = response.content
            
            # Check file size (max 10MB)
            if len(content) > 10 * 1024 * 1024:
                return None
            
            # Save image
            with open(file_path, 'wb') as f:
                f.write(content)
            
            # Validate image
            try:
                with Image.open(file_path) as img:
                    width, height = img.size
                    if width < 50 or height < 50:
                        file_path.unlink()
                        return None
            except Exception:
                file_path.unlink()
                return None
            
            logger.info(f"Downloaded image: {filename}")
            return str(file_path)
        
        return None
    
    except Exception as e:
        logger.error(f"Error downloading image {image_url}: {e}")
        return None


def scrape_website(url: str, output_dir: str = "scraped_data") -> Dict:
    """
    Master function to scrape a website and save all content to a unique directory
    
    Args:
        url: The URL to scrape
        output_dir: Base directory for output (default: "scraped_data")
        
    Returns:
        Dictionary with scraping results and metadata
    """
    
    # Create unique output directory
    unique_dir = create_unique_output_dir(url, output_dir)
    
    result = {
        "url": url,
        "output_directory": str(unique_dir),
        "timestamp": datetime.now().isoformat(),
        "success": False,
        "title": "",
        "images": [],
        "links": [],
        "error": None,
        "files_created": []
    }
    
    try:
        # Setup session
        session = requests.Session()
        session.headers.update({
            "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36"
        })
        
        # Fetch the main page
        response = session.get(url, timeout=30)
        if response.status_code != 200:
            result["error"] = f"HTTP {response.status_code}"
            return result
        
        html_content = response.text
        soup = BeautifulSoup(html_content, 'lxml')
        
        # Extract title
        title_tag = soup.find('title')
        title = title_tag.get_text().strip() if title_tag else "Untitled"
        result["title"] = title
        
        # Save HTML
        html_file = unique_dir / "page.html"
        with open(html_file, 'w', encoding='utf-8') as f:
            f.write(html_content)
        result["files_created"].append(str(html_file))
        
        # Extract and save text content
        # Create a copy for text extraction
        text_soup = BeautifulSoup(html_content, 'lxml')
        for tag in text_soup(['script', 'style', 'nav', 'footer', 'header']):
            tag.decompose()
        
        text_content = text_soup.get_text()
        text_content = re.sub(r'\n\s*\n', '\n\n', text_content).strip()
        
        text_file = unique_dir / "content.txt"
        with open(text_file, 'w', encoding='utf-8') as f:
            f.write(text_content)
        result["files_created"].append(str(text_file))
        
        # Extract and download images
        img_tags = soup.find_all('img')
        downloaded_images = []
        
        for img in img_tags:
            img_src = img.get('src') or img.get('data-src', '')
            if img_src:
                local_path = download_image(img_src, url, unique_dir, session)
                if local_path:
                    image_info = {
                        "src": img_src,
                        "alt": img.get('alt', ''),
                        "local_path": local_path
                    }
                    downloaded_images.append(image_info)
                    result["files_created"].append(local_path)
        
        result["images"] = downloaded_images
        
        # Extract links
        links = []
        for link in soup.find_all('a', href=True):
            link_url = urljoin(url, link['href'])
            link_text = link.get_text().strip()
            links.append({
                "url": link_url,
                "text": link_text
            })
        
        result["links"] = links[:50]  # Limit to first 50 links
        
        # Extract metadata
        metadata = {
            "scraped_at": datetime.now().isoformat(),
            "url": url,
            "title": title,
            "images_count": len(downloaded_images),
            "links_count": len(links),
            "content_length": len(text_content),
            "html_size": len(html_content)
        }
        
        # Add meta tags
        for meta in soup.find_all('meta'):
            name = meta.get('name') or meta.get('property') or meta.get('http-equiv')
            content = meta.get('content')
            if name and content:
                metadata[f"meta_{name}"] = content
        
        # Save metadata
        metadata_file = unique_dir / "metadata.json"
        with open(metadata_file, 'w', encoding='utf-8') as f:
            json.dump(metadata, f, indent=2)
        result["files_created"].append(str(metadata_file))
        
        # Save complete result
        result_file = unique_dir / "scraping_result.json"
        with open(result_file, 'w', encoding='utf-8') as f:
            json.dump(result, f, indent=2, default=str)
        result["files_created"].append(str(result_file))
        
        result["success"] = True
        logger.info(f"Successfully scraped {url} to {unique_dir}")
        
        # Close session
        session.close()
        
    except Exception as e:
        result["error"] = str(e)
        logger.error(f"Error scraping {url}: {e}")
    
    return result


# Example usage and testing
if __name__ == "__main__":
    # Test the scraper
    test_url = "https://example.com"
    result = scrape_website(test_url)
    
    print(f"Scraping result: {result['success']}")
    print(f"Output directory: {result['output_directory']}")
    print(f"Images downloaded: {len(result['images'])}")
    print(f"Files created: {len(result['files_created'])}")
    
    if result['error']:
        print(f"Error: {result['error']}")
    else:
        print("Scraping completed successfully!")
