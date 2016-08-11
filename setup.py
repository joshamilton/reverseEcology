from setuptools import setup, find_packages

def readme():
    with open('README.md') as f:
        return f.read()
    
setup(name='reverseEcology',
      version='0.1',
      description='Reverse ecology analysis of metabolic network reconstructions',
       classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2.7',
      ],
      url='https://github.com/joshamilton/reverseEcology',
      author='Joshua J. Hamilton',
      author_email='joshamilton@gmail.com',
      license='MIT',
      packages=find_packages(),
      include_package_data=True)
