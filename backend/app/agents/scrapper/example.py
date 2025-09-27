#!/usr/bin/env python3
"""
Example usage of the simple web scraper
"""

import sys
from master_scraper import scrape_website


def main():
    """Example usage of the scraper"""
    
    # Check if URL is provided
    if len(sys.argv) < 2:
        print("Usage: python example.py <URL> [output_directory]")
        print("Example: python example.py https://example.com")
        sys.exit(1)
    
    url = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else "scraped_data"
    
    print(f"Scraping: {url}")
    print(f"Output directory: {output_dir}")
    print("-" * 50)
    
    # Scrape the website
    result = scrape_website(url, output_dir)
    
    # Print results
    print(f"\nScraping Result:")
    print(f"Success: {result['success']}")
    print(f"Title: {result['title']}")
    print(f"Output Directory: {result['output_directory']}")
    print(f"Images Downloaded: {len(result['images'])}")
    print(f"Links Found: {len(result['links'])}")
    print(f"Files Created: {len(result['files_created'])}")
    
    if result['error']:
        print(f"Error: {result['error']}")
    else:
        print("\nFiles created:")
        for file_path in result['files_created']:
            print(f"  - {file_path}")
        
        if result['images']:
            print(f"\nImages downloaded:")
            for img in result['images']:
                print(f"  - {img['src']} -> {img['local_path']}")


if __name__ == "__main__":
    main()
