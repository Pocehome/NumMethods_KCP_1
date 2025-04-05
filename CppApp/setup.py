from setuptools import setup, Extension
import pybind11

pybind11_include = pybind11.get_include()

ext_modules = [
    Extension(
        'KSR1_Module',
        ['TaskModule.cpp', 'TestTask.cpp', 'KSR_Task.cpp'],
        include_dirs=[pybind11_include],
        language='c++',
    ),
]

setup(
    name='KSR1_Module',
    version='0.1',
    ext_modules=ext_modules,
    zip_safe=False,
)
