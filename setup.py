from setuptools import setup, find_packages

def readme():
    with open('README.md') as f:
        return f.read()
    
setup(name='reverseEcology',
      version='1.0.0',
      description='Reverse ecology analysis of metabolic network reconstructions',
      classifiers=[
        'Development Status :: 7 - Inactive',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2.7',
      ],
      url='https://github.com/joshamilton/reverseEcology',
      author='Joshua J. Hamilton',
      author_email='joshamilton@gmail.com',
      license='MIT',
      packages=find_packages(),
      install_requires=[
          'cobra>=0.5.4',
          'matplotlib',
          'networkx>=1.11',
          'numpy',
          'pandas',
          ],
      include_package_data=True)
