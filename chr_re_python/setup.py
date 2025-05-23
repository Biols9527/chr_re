from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as f:
    requirements = f.read().splitlines()

setup(
    name="chr_re",
    version="0.1.0",
    author="[Your Name/Organization]",
    author_email="[your_email@example.com]",
    description="A Python framework for phylogenetic character reconstruction and evolutionary analysis.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="[https://github.com/yourusername/chr_re_python]", # Replace with your repo URL
    packages=find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha", # Or "4 - Beta", "5 - Production/Stable"
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License", # Choose your license
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.8",
    install_requires=requirements,
    extras_require={
        "dev": ["pytest>=6.0", "flake8", "black", "mypy"],
        "docs": ["sphinx", "sphinx-rtd-theme"], # Example for documentation
    },
    entry_points={
        "console_scripts": [
            # If you have any command-line runnable scripts
            # "chr_re_cli=chr_re.cli:main",
        ],
    },
    project_urls={
        "Bug Tracker": "[https://github.com/yourusername/chr_re_python/issues]",
        "Documentation": "[https://yourusername.github.io/chr_re_python/]", # If you host docs
        "Source Code": "[https://github.com/yourusername/chr_re_python]",
    },
)
