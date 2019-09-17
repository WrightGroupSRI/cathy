from setuptools import setup, find_packages

setup(
    name="cathy",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "click",
        "click_log",
        "enlighten",
        "catheter_utils @ git+http://panoptes.sri.utoronto.ca:8088/wright-group/catheter_utils.git@a832f6e171ac29e3fae2b6254219b8027c50320e",
        "catheter_ukf @ git+http://panoptes.sri.utoronto.ca:8088/wright-group/catheter_ukf.git@11a1646675d9874ceafc97b5aa93371756619d42",
    ],
    entry_points={"console_scripts": ["cathy=cathy.__main__:main"]},
)
