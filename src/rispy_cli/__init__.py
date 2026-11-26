from .download_data import (
    missing_required_files,
    download_and_extract,
    missing_deepvariant,
    download_and_extract_deepvariant,
)

# Automatically download missing large data on first import
if missing_required_files():
    download_and_extract()

# Automatically download DeepVariant model chunks on first import
if missing_deepvariant():
    download_and_extract_deepvariant()