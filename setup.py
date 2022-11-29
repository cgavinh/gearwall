from setuptools import setup, find_packages

with open('README.md') as readme_file:
    README = readme_file.read()

with open('HISTORY.md') as history_file:
    HISTORY = history_file.read()

setup_args = dict(
    name='gearwall',
    version='0.1.0',
    description='My personal code, including simulation software, data analysis, and general utilities',
    long_description_content_type="text/markdown",
    long_description=README + '\n\n' + HISTORY,
    license='MIT',
    # Which Python importable modules should be included when your package is installed
    # Handled automatically by setuptools. Use 'exclude' to prevent some specific
    # subpackage(s) from being added, if needed
    packages=find_packages(),
    author='Coire Gavin-Hanner',
    author_email='c.gavin.hanner@gmail.com',
    keywords=['Data', 'Analysis', 'Quantum'],
    url='https://github.com/cgavinh/gearwall',
    download_url='https://pypi.org/project/gearwall/',
    include_package_data=True
)


if __name__ == '__main__':
    setup(**setup_args)