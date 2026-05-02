default:
    @just --list

# Open nb01 in a live marimo kernel with deps provisioned from its PEP 723 header.
notebook:
    @[ -f .env ] || (echo "Missing .env. Copy .env.example and paste your FINNGENIE_TOKEN." && exit 1)
    env -u PYTHONPATH uvx marimo edit --sandbox notebooks/nb01_pcsk9_walkthrough.py
