# Simple Web Scraper

A simple, focused web scraper that takes a single URL and downloads all content (HTML, text, images) to a unique output directory.

## Features

- **Single Function Interface**: One master function `scrape_website(url, output_dir)`
- **Unique Output Directories**: Each scraping session gets its own timestamped directory
- **Complete Content Download**: HTML, extracted text, images, and metadata
- **Image Validation**: Downloads and validates images (size, format, dimensions)
- **Link Extraction**: Extracts all links from the page
- **Metadata Collection**: Saves page metadata including meta tags
- **Error Handling**: Comprehensive error handling and logging

## Installation

Install the required dependencies from the backend directory:

```bash
# From the backend directory
pip install -r requirements.txt

# Or with virtual environment
source venv/bin/activate
pip install -r requirements.txt
```

## Usage

### As a Python Module

```python
from master_scraper import scrape_website

# Scrape a website
result = scrape_website("https://example.com", "my_output_folder")

# Check results
if result['success']:
    print(f"Scraped successfully to: {result['output_directory']}")
    print(f"Downloaded {len(result['images'])} images")
    print(f"Found {len(result['links'])} links")
else:
    print(f"Error: {result['error']}")
```

### Command Line Usage

```bash
# Basic usage
python example.py https://example.com

# With custom output directory
python example.py https://example.com my_custom_output

# Example with a real website
python example.py https://httpbin.org/html
```

## Function Signature

```python
def scrape_website(url: str, output_dir: str = "scraped_data") -> Dict:
    """
    Master function to scrape a website and save all content to a unique directory
    
    Args:
        url: The URL to scrape
        output_dir: Base directory for output (default: "scraped_data")
        
    Returns:
        Dictionary with scraping results and metadata
    """
```

## Return Value

The function returns a dictionary with the following structure:

```python
{
    "url": "https://example.com",
    "output_directory": "/path/to/unique/directory",
    "timestamp": "2025-08-17T12:34:56.789012",
    "success": True,
    "title": "Page Title",
    "images": [
        {
            "src": "image_url",
            "alt": "alt text",
            "local_path": "/path/to/downloaded/image"
        }
    ],
    "links": [
        {
            "url": "link_url",
            "text": "link text"
        }
    ],
    "error": None,
    "files_created": [
        "/path/to/page.html",
        "/path/to/content.txt",
        "/path/to/metadata.json",
        "/path/to/images/image1.jpg"
    ]
}
```

## Output Directory Structure

Each scraping session creates a unique directory with the following structure:

```
scraped_data/
└── example.com_20250817_123456_a1b2c3d4/
    ├── page.html              # Original HTML
    ├── content.txt           # Extracted text content
    ├── metadata.json         # Page metadata
    ├── scraping_result.json  # Complete scraping results
    └── images/               # Downloaded images
        ├── image1.jpg
        ├── image2.png
        └── ...
```

## Configuration

The scraper has built-in sensible defaults:

- **Image Size Limit**: 10MB per image
- **Minimum Image Dimensions**: 50x50 pixels
- **Supported Image Formats**: .jpg, .jpeg, .png, .gif, .webp, .svg
- **Request Timeout**: 30 seconds
- **Image Download Timeout**: 10 seconds per image

## Error Handling

The function includes comprehensive error handling:

- Network timeouts
- Invalid URLs
- HTTP error codes
- Image validation failures
- File system errors

All errors are logged and returned in the result dictionary.

## Use as a Tool

This scraper is designed to be used as a tool in larger systems. The single function interface makes it easy to integrate:

```python
# Use in an agent or automation system
def scrape_tool(url: str) -> str:
    result = scrape_website(url)
    if result['success']:
        return f"Successfully scraped {url} to {result['output_directory']}"
    else:
        return f"Failed to scrape {url}: {result['error']}"
```

## Dependencies

- `requests`: HTTP requests
- `beautifulsoup4`: HTML parsing
- `Pillow`: Image validation
- `lxml`: Fast XML/HTML parsing

## Example Output

```
Scraping: https://example.com
Output directory: scraped_data
--------------------------------------------------

Scraping Result:
Success: True
Title: Example Domain
Output Directory: scraped_data/example.com_20250817_123456_a1b2c3d4
Images Downloaded: 0
Links Found: 1
Files Created: 4

Files created:
  - scraped_data/example.com_20250817_123456_a1b2c3d4/page.html
  - scraped_data/example.com_20250817_123456_a1b2c3d4/content.txt
  - scraped_data/example.com_20250817_123456_a1b2c3d4/metadata.json
  - scraped_data/example.com_20250817_123456_a1b2c3d4/scraping_result.json
```
