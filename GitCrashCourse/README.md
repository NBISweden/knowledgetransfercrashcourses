
# A crash course to `git`

`git` is a version control system that can  be used to keep track
of changes made to files during the development of a project. The
version-controlled files are stored in a central _repository_. The
repository can be shared with others, either as a _public_ repository
or as a with selected users as a _private_ repository. Several users
can participate in the development simultaneously and independently.

`git` projects are often the development of code of some software,
but can also be, e.g., manuscript writing, typically in LaTex
(but maybe not in, e.g., MS word, see next subsection).

__Files for which `git` version control does not work so well__

- _MS Word, Excel, MacOSX Pages and Numbers, pdf etc._ `git` is really
good at identifying individual changes (e.g., a word or a line)
between old and new pure text files and storing these efficiently.
This does not work when the files are not pure text, typically
meaning that the whole file is saved at each change. However,
including these type of files meant as documentation only (e.g.,
platform reports) and not being subject to frequent changes of
course works.
- _Very large files, typically data files._ One should avoid
adding very large files to a git repo as this might crash the
repo. Formally, files over 50Mb are considered very large;
however, it's good practice to try to keep file sizes much
smaller. There are ways around this if really needed, but
large files should generally be avoided.

I will below assume that you work on a computer where you can
open a unix-like terminal (either your laptop or UPPMAX). In
fact, all commands in the rest of this text assumes that you
work in such a terminal and that you have _some_ experience
using it.

## What about Github and Bitbucket?
The central repository (or _repo_ for short) can be on a local
server or laptop, but it is always much better to keep it at
some repo service provider. There are several such providers,
Github, Bitbucket, Gitlab, etc.; the choice of which to use is
almost arbitrary (google to find out more).
_NBIS conflicts_ has an professional account with Bitbucket,
which is why we use it in this project.

You will, most likely, need to apply for an an account to the
central repository you choose to use. This is free and will
provide you with your own user-name and a password.

## Installing `git`

### UPPMAX
`git` version 1.8.3.1 is installed on rackham; I usually use
this version, but it is possible to get newer versions using
the `module` system. Type
```
module spider git
```
to see what versions are available, in this case.

### MacOSX
I think git is installed by default in MacOSX. Possibly
one need to install Xcode and its command line tools.

Otherwise, see
[here](https://www.atlassian.com/git/tutorials/install-git) for
a tutorial on `git` installation (scroll down to MacOSX).

### Windows
I don't have much experience of Windows, but if you have an
application providing unix terminal, maybe this also provides ’git’.

Otherwise, see
[here](https://www.atlassian.com/git/tutorials/install-git) for
a tutorial on `git` installation (scroll down to Windows).

### GUI version?
This crash course is for working with `git` in the terminal,
but there are also graphical user interfaces (GUIs) for git.
I have never used these and don't really know anything about
this (google if you are interested).

## Most important git commands
### Cloning a repo `git clone`
We will consider the situation where you want to clone an existing
repository and create a working directory (_wd_) on your computer.
Typically you have a web-address to the repository; we'll assume
here that it is a Bitbucket repo, but it works similarly for
Github, etc.

1. open a unix terminal and `cd` to the place where you want
place the _wd_.
2. Open the address in you browser and click the button
`clone` in the upper right corner.
3. In the window that opens, click the button that says `SSH`
in the top right corner and change to `HTTPS` (unless it
  already says `HTTPS` there).
4.Now copy the
`git clone https://yourusername@bitbucket.org/scilifelab-lts/reponame.git`
text and paste it into your terminal, i.e., it should say
```
git clone https://yourusername@bitbucket.org/scilifelab-lts/reponame.git
```
in your terminal. Click enter.

This creates a folder with the same name as the repo and downloads
all files in the repo into that folder. It also sets up the folder
as a `git` working directory with a _local_ `git` project. The
Bitbucket repo is referred to as a _remote_.

### Updating a repo `git pull`
The repo might be changed after you have cloned it. To get all
new changes to your _wd_ you must explicitly update it by typing
```
git pull
```
(you must be somewhere inside your _wd_ when doing this).
This will probably require that you type your Bitbucket password.

This will update all files that are out of date in your _wd_.
If you have done some changes to your files, `git` will usually
warn about this and not do the pull. It is good to to take care
of such changes (see below how) before `git pull`.

### Do I have any changes? `git status` and `git diff`
To check for changes type
```
git status
```
This will write a small report, including
- if you have any _Changes not staged for commit_,i.e., changed
files that have not been _committed_ to `git` yet.
- if you have any _Untracked files_, i.e., new files that git
doesn't know about.

For how to handle these, see below on `git add` and `git commit`.

If you want to see what the _uncommmitted_ changes are for a file
(maybe so you can decide whether it should be kept or not), you can
type
```
git diff path/to/files
```
This will show line-by-line differences between the file and the last committed version of it.

### Tell `git` about your changes: `git commit` and `git add`
If you have made changes to a file that you want to keep in
the `git` repo, you should _commit_ them to `git`:
```
git commit -m "brief description of changes" path/to/files
```
This will include the changes in the local git project.

If you have a new file that you want to add to `git`, you
can type
```
git add path/to/file
```
You then need to also commit this file as described above
to complete the file addition.

### How to remove changes or files: `git checkout`, `git mv` and `git rm`
If you have changes that you don't want to keep you can
remove them from the file manually or just checkout the latest
version of the file in git, by typing
```
git checkout path/to/file
```
This will remove all changes you have made in the file.

If you want to change the name of a file in git or just move it,
you should always use `git mv` instead of the system command `mv`
```
git mv path/to/oldfilename path/to/newfilename
```
If you have a whole file that is no longer needed and you
want to remove it from git, similarly use `git rm` and not
the system `rm`:
```
git rm path/to/file
```

Finally, if you want to remove an _untracked file_, use system `rm`
```
rm path/to/file
```

### Letting the _remote_ know: `git push`
To update the _remote_ repo with your committed changes, you need
to _push_ to the _remote_
```
git _push
```
This will probably require that you type your Bitbucket password.

#### Merge conflicts
When several developers are working on the same repo, they might
add changes to the same file (almost) at the same time. `git` will
try to resolve them and merge the two version automatically -- you
might be asked to cofirm the merging (usually a text editor window
pops up; you usually just have to close it).
Sometimes, however, there might be conflicts between your changes
and someone else's changes. In that case you need to resolve these
conflicts and merge them before being able to push your changes.
When `git` reports that there are merge conflicts for a file,
you should:
1. Open the file in a text editor
2. Find the merge conflicts in the file; these are marked (by
`git`) using sequences of ``>>>>>>``, ``<<<<<<<<``, and ``=======`` --
usually just a few lines differ.
3. Try to figure out if you can manually merge the two version --
maybe you have tried to fix the same thing, then decide which
solution is best or maybe the conflict is really minor.
4. Fix it, remove the conflict markers and the text that is
not needed (remember to check if there are further conflicts
in the file)
5. Close the file, `git commit` it and `git push`
6. If the conflicts are hard to resolve, maybe it is unclear what
you can remove or maybe very large parts of the file is affected,
don't hesitate to ask for help -- you might need to communicate
with the person who committed the text in conflict with yours or
you might want to consult someone more experienced.


## Create a new branch
If you want to introduce modifications or develop new features for a certain application without having to change anything in the main branch `master`,
then you can create a new branch.


Let's say you are on the `master` branch and want to create a new branch. Here's what you'll do:
```
git checkout -b <my branch name>.
```

This command will automatically create a new branch and then move you to that branch.
To confirm that your new branch was created you can use the `git branch` command.
The branch name with the asterisk next to it indicates which branch you're pointed to at that given time.


Once you've created a branch, you can save your progress with `git add` and `git commit` as usual.
Remember that executing `git push` will upload the changes _remotely_ to the current branch `<my branch name>`
(_NB!_ You may perhaps get a message from git and  be asked to push explicitly to the remote version of the branch -- I usually just follow the syntax that is suggested in that message).


Now if you are ready to make all the modifications in the branch to be part of the main branch, then you can merge it
to `master` by executing:

```
git checkout master
git merge <my branch name>
```

Note that at any point in time, you can switch to the `master` branch by `git checkout master`.
Similary, if you want to switch to a different branch then do `git checkout <my branch name>`.




## Good `git` practices
When developing code/text (i.e., changing files) in `git`, there
are some things one can do to avoid _merge conflicts_:
- Start each developing session by `git pull` so you don't work
on outdated files.
- End each session by `git push` for much the same reason.
- If you feel you need to do mayor changes that may take longer
time and maybe lots of try-outs, you should read up on
_branching in git_ and/or _forking repositories_ (above, or google it).

## Further reading
I tried to make this crash course rather NBIS project-specific.
There are several other, more general, crash courses or tutorials
that complements, improves or extends this text. You can find
them by googling.

Github provides a rather good `git` documentation which can be googled
Specifically, one can look at the commands `git log`, `git stash`,  
_branching in git_, and the concepts of _forking repositories_ --
but these belong to the advanced course:)
