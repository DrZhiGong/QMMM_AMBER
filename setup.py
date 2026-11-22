"""Tell pip how to install the qmmm command and the files it uses."""

from setuptools import setup

setup(
    name="amber-orca-qmmm",
    version="0.1.0",
    description="Relax an AMBER environment around a fixed ORCA QM region",
    python_requires=">=3.11",
    packages=["AMBER_ORCA", "AMBER_ORCA.src"],
    package_dir={"AMBER_ORCA": "."},
    entry_points={"console_scripts": ["qmmm = AMBER_ORCA.qmmm:cli"]},
    package_data={
        "AMBER_ORCA": [
            "README.md", "settings.json", "run.sh", "src/runtime/*.sh", "dat/*", "doc/*.md",
            "examples/*/README.md", "examples/*/*.json", "examples/*/*.prmtop",
            "examples/*/*.rst7", "examples/*/*.in",
        ],
    },
)
