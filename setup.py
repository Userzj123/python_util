from setuptools import setup
import glob

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="python_utils",
    version="0.1.0",
    packages=["turb", 'pyutils'],
    package_dir = {'': 'src'},
    py_modules = ['turb.channel_plot', 'pyutils.plot_utils', 'pyutils.data_utils'],
    # scripts=glob.glob('src/turb/*.py'),
    install_requires=[
        'Pillow>=9.5.0',
        'imageio>=2.31.1',
        'numpy>=1.25.0',
        'matplotlib>=3.7.1'
        ],
    license="MIT",
    url="https://github.com/Userzj123/python_util.git",
    author="Zejian You",
    author_email="userzj123@gmail.com",
    description="Utilities for python in research.",
    long_description=long_description,
    long_description_content_type="text/x-rst",
    test_suite="nose.collector",
    tests_require=["nose"],
)