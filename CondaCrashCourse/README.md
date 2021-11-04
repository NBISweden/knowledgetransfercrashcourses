# Conda Crash Course

**Todo:**

- It would be neat if one could create a conda environment for a decent,
*intuitive* cross-platform text-editor that have capabilities (correct
indentation, code-formatting, etc.) to properly handle yaml- and
snakemake-formatted files, at least, but preferably also markdown and maybe
also python and R/Rmarkdown. `Atom` and `Sublime` would be alternatives, but
there are no conda-packages for these yet (Sublime is formally 'licensed', but atom is open-source ; MIT-license). Maybe a TODO could be to create
an `atom` or `sublime` package
-- or find an alternative editor. (`Emacs` is actually my editor of preference,
but that fails the *intuitive* criterium... and so does ’vim’)
- need feedback on how this works with Windows (i.e., Windows unix-emulator
or Windows 10)
- We could add info on checksums, but I think this is an overkill?


## Introduction

Conda is a cross-platform system for software package management.
It greatly simplifies software installation on your computer.
Moreover, it allows, in a very simple way, creating project-specific
environments with a selection of installed programs needed for a particular
project. Each time work is performed on the project, the environent is
activated, and when work is done done, the environment is deactivated the
standard software environment for system is back.

Conda works on MacOS, Windows and Unix/Linux systems. This is handy for
reproducibility as, in principle, a script developed inside a conda environment
on, say, a Mac laptop can be transferred to, say, a linux system on UPPMAX, the
corresponding conda environment for linux can be created and the script can be
used there as well.

Conda packages of programs are available from different *channels*, including
[anaconda](https://repo.anaconda.com) (often set as the `defaults` or `main`
channel in conda), [conda-forge](https://conda-forge.org) and
[bioconda](https://bioconda.github.io); all (or at least the majority of)
channels are isted in the [anaconda cloud](https://anaconda.org)). A large part
of the available open-source bioinformatics programs are covered by these
channels. Commercial or other programs  with restricted license can typically
not be provided by these channels -- however, for reproducibility, open source
programs are a more natural choice.

It should be mentioned that, while `conda` is available for MacOS, Windows and
Unix/Linux, certain packages may not be available for all platforms. This could
be due to the actual program being platform-specific, to difficulties to port
the program between platforms or to channel policies (*bioconda* not supporting
Windows?).

## Additional info
There's also nice information and very good tutorials about conda available in the
**Tools in reproducible research** [course material](https://nbis-reproducible-research.readthedocs.io/en/latest/conda/).


## Installation

Conda can be installed using the Miniconda installer. The installer is
downloaded from [here](https://docs.conda.io/en/latest/miniconda.html); choose
the installer appropriate for your operating system and download it.

Detailed instruction can be found [here](https://conda.io/projects/conda/en/latest/user-guide/install/index.html),
but briefly, the downloaded file is a *bash* script or a Windows executable,
respectably, which is executed and the user is guided through the installation.

Once you've followed the installation instructions the base conda environment
will be activated by default when you open up a new terminal window. You can
see this by looking at your prompt which should be prefixed with `(base)` or
`(anaconda3)`. You can deactivate the base environment by running
`conda deactivate` and re-activate it by running `conda activate`. Notice how
your prompt changes as you do this. See [more below](#activating-and-deactivating-an-environment)
about activating and deactivating environments.

### Updating conda

Conda can be used to update itself, simply by typing:

```
conda update conda
```

### Local UPPMAX conda

On the UPPMAX clusters (e.g., )*rackham*, *bianca*, *snowy* and *irma*),
`conda` is also available through the UPPMAX module system. This module provides
a *local* variant of `conda`, that i, it accesses local mirrors of the major
bioinformatics-relevant *conda* channels (including *anaconda*, *bioconda*, and
*conda-forge*; for a full listing do `conda config show channels`).

If your project involves human sensitive data and, hence, work on the secure
cluster *bianca*, you will need to use this local `conda` module, as *bianca*
does not have any internet connections. The UPPMAX local *conda* can be loaded
(while on an UPPMAX cluster) simply by typing:

`module load conda`

Unloading the module is done by:

`module unload conda`

## Usage

(Note! If conda seems to be extremely slow, read [this section first](#mamba----if-your-conda-is-very-slow-and-shaky-mamba))

### Finding a conda packages

To find out if there exists a conda packages for a certain package, the easiest
way is actually to google `conda <package name>`. You can also search
[Anaconda cloud](https://anaconda.org/) for packages
directly. Pages for packages at *Anaconda cloud* tells:

- which channel provides the packages
- which systems the package is provided for, and
- which version of the software the package provides (other versions can be
found under the *files* tab or as described below).

`conda` itself can also be used to search for packages:

```
conda search --channel <channel> <package>
```

A drawback is that you need to tell conda what channel to search using the
option `--channel` (although the option can be repeated in the same command to
search in several channels at the same time, e.g.,
`conda search --channel bioconda --channel conda-forge <package>`).

#### Searching on UPPMAX

On UPPMAX, `conda` is set up to search *all* local channel mirrors, by simply
typing:

```
conda search <package>
```

### Creating environment

Conda environments can be created indifferent ways. Since reproducibility is a
major theme for these crash courses, we will here focus on the approach where a
conda environment file is used. Separate conda environments files are used for
separate environments. The environment file should be in `yaml` format and can
look like this:

*file: smthg.yaml*
```
Name: smthg
# indentations are important and should be made up of spaces (soft-tabs).
channels:
  - anaconda
  - conda-forge
dependencies:
  - smthg =1.0.1
  - smthelse =1.5.3
```

Breaking down the content:
- The `name` is the name you want to call the conda environment. This could be,
e.g., the name of the contained software or of the project it will be used in.
- Under `channels`, the channels that supply the required conda packages are
listed. Several channels can be given, but it is good practice to only list
those necessary -- avoid unnecessary safe-guarding!
- Finally, the required conda packages are listed under `dependencies`. Also
here, several packages can be listed in the same file. Conda usually makes a
good job in taking care of any incompatibilities between packages, and warns if
it fails doing this.   
For reproducibility reasons it is good to also list the version of the software
required. This is done by adding a ’=’ sign followed by the version (to find
what versions are  as packages, `conda search` can be used, see above).

The conda environment is then created by

```
conda env create --file smthg.yaml
```

This will ask you to verify the installation, then perform it and report
success or failure. Reasons for failure is most often due to:

- the name of the environment is already in use -- change environment name or
remove the old environment (see below).
- A package is not available from the channel or not available for the
platform (i.e., MacOS, Windows or Linux/Unix) -- correct channel or use another
platform (e.g., UPPMAX).
- Versions of two or more packages, as required in the conda environment yaml
file, are not compatible -- try with other versions. For example, a common
incompatibility is that one package is built for `python 2` and another is built
 for `python 3` -- these are not compatible.

When setting up an environment for a new project it can be a good idea to
**not** set explicit versions for packages, unless you know beforehand that your
project requires a specific version. This way you can let conda attempt to
install the most up to date compatible versions for the listed packages. Once
the environment is created you can update the environment file with the exact
versions that were installed (use `conda list` to see package versions for the
activated environment).

##### A note on editing yaml-files
A really annoying drawback with the yaml-format is that indentations
_must_ be as spaces and not as ASCII tab-characters_. However, many
(but not all) text editors handles this `tabs-as-spaces` requirement
quite well; that is, when you hit the tab-key it will insert a sequence
of spaces instead of the ASCII tab character
([Atom](https://github.com/atom/atom]) and
[Sublime](https://www.sublimetext.com) are quite easy, "free",
multi-platform editors that can be set to handle this issue).

### Activating and deactivating an environment

A succesfully installed conda environment is activated by typing:

```
conda activate smthg
```

However, often `conda` will complain that it is not initiated for bash. this is
done by typing:

```
conda init bash
```

Then retry the activation command. Now you will have access to the programs
provided by the conda environment. You will also have access to the standard
programs of your system, but if you have the same program in your standard
system and in your conda environment, the one in the conda environment will be
used.

The conda environment is deactivated by typing:

```
conda deactivate
```

The standard environment is then restored and the conda environment's software
is not accessible.

### Updating an environment

To add additional packages to an environment:

1. First update the conda environment files in question by adding additional
packages (and channels as needed) or change the version of existing packages.
2. Then type:

```
conda env update --file smthg.yaml
```

and the environment will be updated.


### Deleting an environment

Sometimes an outdated environment needs to be deleted. This is done by typing:

```
conda env remove --name smthg
```

`conda` then asks for verification, then deletes the environment and reports
back.

### Keeping your conda clean

By default, conda saves all downloaded software tarballs (archives of source
code), redundant packages (e.g., older, replaced versions) and index-cache
because, in case they are needed, it might save a little time not having to
download/create them again. However, after some time all these saved files will
make the conda-installation cluttered, taking up disk-space and possibly
affecting performance. It is there a good practice to clean the conda
installation from time to time. This done by typing:

```
conda clean --all
```

`conda` asks for verification for cleaning, then perform the cleaning and
reports back.

### Mamba -- if your conda is very slow and shaky

As the package channels has grown very large, it has turned out that the
original conda implementation, at times, has problems parsing all packages.
It can therefore be very slow and maybe have problems resolving package
conflicts properly. This issue is being worked on and may, by the time you
read this, already be solved.

However, a very good solution for the time being is to use `mamba`, which
itself is a conda -package and which solves the problems mentioned above.
You install `mamba` using conda; we recommend that you install it in your
`base` environment. To do this, ascertain that you have your `base` conda
environment (and no other environment) open and type:

```
conda install mamba
```

If you prefer, you can install `mamba` in a dedicated environment, as
described above; just remember to activate this environment every time
you want install packages using `mamba`.

You then simply substitute `mamba` for `conda` in almost all conda commands
described above (adn below in the exercises), e.g.,

```
mamba env create -f smthg.yaml
```

Some exceptions are `conda activate`, `conda deactivate` and `conda init`.

## Exercises

1. Create a conda environment for snakemake using the file [snakemake.yaml](snakemake.yaml).
Activate the environment and test that you have access to the program by typing:

```
snakemake snakemake --version
```

Then deactivate the environment.  
Did it work? Did you get the right version?

2. Use a text-editor and create your own conda environment file `bedtools.yaml`, google or conda search to find a package for bedtools
(you can select some version x.y.z. if you like) and add that to the
environment file.
Create the environmentand activate it. Check that you can run bedtools.
3. Add `samtools` to the environment file and update the environment. Check
that it worked.
4. Deactivate and remove the environment.
