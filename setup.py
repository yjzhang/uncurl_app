from setuptools import setup, find_packages

# TODO: this won't work...
install_requires = []

setup(
    name='uncurl_app',
    version='0.3.0',
    author='Yue Zhang',
    author_email='yjzhang@cs.washington.edu',
    url='https://github.com/yjzhang/uncurl_app',
    license='MIT',
    scripts = [
            'uncurl_app_split_seq'
    ],
    install_requires=install_requires,
    packages=find_packages("."),
    include_package_data=True,
    test_suite='nose.collector',
    tests_require=['nose', 'flaky'],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.6',
    ],

)

