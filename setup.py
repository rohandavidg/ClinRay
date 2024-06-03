from setuptools import setup
import os
with open(os.path.join(os.path.dirname(__file__), 'bin', 'VERSION')) as version_file:
    version = version_file.read().strip()

with open('requirements.txt') as requires_file:
    requires = requires_file.read().split('\n')

setup(
    name='ClinRay',
    version=version,
    description='ClinRay Homology detection tool',
    url='https://github.com/rohandavidg/ClinRay',
    author='Rohan Gnanaolivu',
    author_email='gnanaolivu.rohandavidg@mayo.edu',
    license='MIT',
    packages=['bed', 'bam'],
    scripts=[
        'bin/prerocessing.py',
        'bin/train.py',
        'bin/model.py',
        'bin/inference.py'
    ],
    install_requires=requires,
    zip_safe=False,
    include_package_data=True
)
