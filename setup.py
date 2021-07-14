from setuptools import setup, find_packages

setup(
    name="cathy",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "click",
        "click_log",
        "enlighten",
        "matplotlib",
        "catheter_utils",
        "catheter_ukf",
        "dicom_utils",
        "dicom_art",
        "get_gt"
    ],
    entry_points={"console_scripts": ["cathy=cathy.__main__:main"]},
)
