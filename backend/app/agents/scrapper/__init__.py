"""
Simple Web Scraper Package

Main function: scrape_website(url, output_dir)
"""

from .master_scraper import scrape_website

__version__ = "1.0.0"
__all__ = ["scrape_website"]
