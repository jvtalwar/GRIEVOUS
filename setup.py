from setuptools import setup, find_packages

setup(
    name = "grievous",
    version = "0.1.0",
    author = "James V. Talwar",
    author_email = "jtalwar@ucsd.edu",
    description = "Generalized Realignment of Innocuous and Essential Variants Otherwise Utilized as Skewed",
    license = "LICENSE",
    url = "https://github.com/jvtalwar/GRIEVOUS", 
    packages = find_packages(),
    python_requires = '>=3.7',
    install_requires = ["pandas >= 1.3.4"],

    scripts=["grievous/generalGrievous.py"],

    entry_points = {
        "console_scripts": [
            "grievous = grievous.generalGrievous:main"
        ]
    }

)