# Git Crash Course

## Why should I use it?
`git` is one of the essential corner-stone for reproducible research.
Not only does it provide a backup for all your code, but it also keeps
track of the changes you make. So in the not uncommon case you make a
lot of inspired changes to a script that in hind-sight turned out not
to be so brilliant, you can easily go back and retrieve your latest
working version of that script and continue from there.

But there is more. You can develop your code on your laptop, and then
when you need to run it large-scale on a cluster, `git` allows you
create a synchronized clone on the cluster and you can continue
developing the code there... or in both places. In fact, if you want
to get advanced, `git` provides means for several collaborators
to develop on the same code on different computers.

Finally, it also mediates sharing the final code publicly , e.g., to
article reviewers etc. It is even possible to automatically create a
nice web-interface to the code.

This course covers only the basic features of `git`. To learn more
advanced issues, toyou can look at some of the

## Introduction

`git` is a version control system that can  be used to keep track
of changes made to files during the development of a project. The
version-controlled files are stored in a central _repository_. The
repository can be shared with others, either as a _public_ repository
or as a with selected users as a _private_ repository. Several users
can participate in the development simultaneously and independently.

`git` projects are often used in the development of code of some
software, but can also be, e.g., manuscript writing, typically in
LaTex (but maybe not in, e.g., MS word, see next subsection).

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
adding very large files to a git repository as this might crash the
repository. Formally, files over 50Mb are considered very large;
however, it's good practice to try to keep file sizes much
smaller. For this reason, data files are typically _not_ included in
git repostiories. There are ways around this if really needed, but
large files should generally be avoided.

I will below assume that you work on a computer where you can
open a unix-like terminal (either your laptop or an UPPMAX login
server). In fact, all commands in the rest of this text assumes that
you work in such a terminal and that you have _some_ experience
using it.

#### Additional info
There's also a very good tutorial about conda available in the
**Tools in reproducible research** [course material](https://uppsala.instructure.com/courses/51980).


## What about GitHub, GitLab and Bitbucket?
The central repository (or _repo_ for short) can be on a local server
or laptop, but it is always much better to keep it at some repository
service provider. There are several such providers, the most common are
GitHub, Bitbucket, and GitLab. All of these offer free accounts that provide
you with your own user-name and a password and allows you to create
repositories or download from public repositories. There are some, more
or less, subtle differences between these providers (see, e.g., discussions
[here](https://www.gangboard.com/blog/github-vs-gitlab-vs-bitbucket)),
but the choice of which to use is almost arbitrary. I will here mainly refer
to the GitHub web interface, but will try to mention differences in Bitbucket.
Unfortunately, I don't have experience of GitLab.

For this crash course, you will not (I believe) need a GitHub account, but
for your future work, you will (almost certainly) need to apply for an an
account to the central repository you choose to use; here follows links to the
account application pages for [GitHub](https://github.com/join),
[Bitbucket](https://id.atlassian.com/signup?application=bitbucket&continue=https://bitbucket.org/account/signin/?optintocst=1&next=/?aidsignup=1),
[GitLab](https://gitlab.com/users/sign_up?test=capabilities).

#### What does a repo look like on its webpage?
The exact design will vary a little among different providers, but typically
you will start at the _code_ or _source_ tab of the repo.

1. At or near the top the repo's _name_ is given.
2. Either below that or in the sidebar, there will be links to
different _repo tabs_ -- you will most often need the _code_/_source_
tab.
3. Below that there is a roll-down menu with the available
[branches]()#create-a-new-branch)in the repo; typically you are interested
in the _master_ branch.
4. Next, the file and folder content of the repository is listed.
5. Below that, the content of the README in the top repo folder is shown.
This README file (typically in _markdown_-format, _.md_) contains a
description of the repo provided by the owner.
6. Lastly, but perhaps most importantly, you will find, somewhere at
the top wright of the page, a roll-down menu button named _Code_ or
_Clone_ that is the key to your access to the repo code (see further
[below](#cloning-a-repository-git-clone))

By clicking on a file (in _5_), you can display its content. Similarly,
clicking on a subfolder opens a new webpage displaying the subfolder's
content essentially in the same manner as described above.

## Installing `git`

### UPPMAX

`git` version 1.8.3.1 is installed on rackham; I usually use
this version, but it is possible to get newer versions using
the `module` system. Type
```
module spider git
```
to see what versions are available, in that case.

### Other Unix/Linux

Most Unix/Linux systems probably provides `git` through their installer
(`apt-get`, `yum`, etc.).
Alternatively, [install using `conda`](#install-using-conda) or see
[here](https://www.atlassian.com/git/tutorials/install-git#linux) for an
installation tutorial (if needed, scroll down to Linux).

### MacOSX
An Apple version of `git` is provided by Apple's Xcode program in its command
line tools. Install Xcode from AppStore; open a terminal and type

```
xcode-select --install
```

For a standard, i.e., non-Apple version, `git`,
[install using `conda`](#install-using-conda) or see
[here](https://www.atlassian.com/git/tutorials/install-git#mac-os-x) for
a tutorial on `git` installation (if needed, scroll down to MacOSX).

### Windows

I don't have much experience of Windows, but if you have an
application providing unix terminal, maybe this also provides ’git’.

Otherwise, [install using `conda`](#install-using-conda) or see
[here](https://www.atlassian.com/git/tutorials/install-git#windows) for
a tutorial on `git` installation (if needed, scroll down to Windows).

### Install using `conda`

Finally, you can install `git` on your system, using `conda` (see the
[CondaCrashCourse exercises](../CondaCrashCourse/README.md#exercises)).

### GUI version?
This crash course is written for working with `git` in the terminal,
but there are also graphical user interfaces (GUIs) for git.
These are not covered by this crash course (`google` if you are interested).

## Most important git commands

### Cloning a repository `git clone`

We will consider the situation where you want to clone an existing
repository and create a git working directory (_gwd_) on your computer.
Typically you have a web-address to the repository; we'll assume
here that it is a Github repository, but it works similarly for
Bitbucket, etc.

1. open a unix terminal and `cd` to the place where you want
place the _gwd_.
2. Open the address in you browser and click the button
`clone` in the upper right corner.
3. In the window that opens, click the button that says `SSH`
in the top right corner and change to `HTTPS` (unless it
  already says `HTTPS` there).
4. Now copy the address and type
`https://github.com/<owner>/<reponame>.git`
text and paste it into your terminal after `git close`, i.e., it should say
```
git clone https://github.com/<owner>/<reponame>.git
```
in your terminal. Click enter.

This creates a folder with the same name as the repository and downloads
all files in the repository into that folder. It also sets up the folder
as a `git` working directory with a _local_ `git` project. The
Bitbucket repository is referred to as a _remote_.

### Updating a repository `git pull`

The content of the remote repository might be changed after you have
cloned it.  Your local _gwd_ does not automatically keep up-to-date
with these changes; instead you must explicitly update your _gwd_
by typing
```
git pull
```
(you must be somewhere inside your _gwd_ when doing this).
This will probably require that you type your account password.

This will update all files that are out of date in your _gwd_.
If you have done some changes to your files, `git` will usually
warn about this and not do the pull. It is good to to take care
of such changes (see below how) before `git pull`.

### Have I made any changes? `git status` and `git diff`
To check for changes type
```
git status
```
This will write a small report, including
- if you have any _Changes not staged for commit_, i.e., changed
files that have not been _committed_ to `git` yet.
- if you have any _Untracked files_, i.e., new files that git
doesn't know about.

For how to handle these, see below on `git add` and `git commit`.

If you want to see what the _uncommmitted_ changes are for a file
(maybe so you can decide whether it should be kept or not), you can
type
```
git diff path/to/file
```
This will show line-by-line differences between the file and the
last committed version of it.

### Tell `git` about your changes: `git commit` and `git add`
If you have made changes to a file that you want to keep in
the `git` repository, you should _commit_ them to `git`:
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

### How to reset, remove, move or rename files: `git checkout`, `git mv` and `git rm`
If you have changes that you don't want to keep you can
remove them from the file manually or just checkout the latest
version of the file in git, by typing
```
git checkout path/to/file
```
This will remove all changes you have made in the file.

If you want to change the name of a file in git or just move it to
another folder in _gwd_, you should **always** use `git mv` instead of
the system command `mv`
```
git mv path/to/oldfilename path/to/newfilename
```
If you have a whole file that is no longer needed and you
want to remove it from git, similarly **always** use `git rm` and not
the system `rm`:
```
git rm path/to/file
```

Finally, if you want to remove an _untracked file_, use system `rm`
```
rm path/to/file
```

### Letting the _remote_ know: `git push`
To update the _remote_ repository with your committed changes, you need
to _push_ to the _remote_
```
git push
```
This will probably require that you type your git account password.

#### Merge conflicts

When several developers are working on the same repository, they might
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
usually just a few lines differ. I usually search for "<<<" in my
text editor to locate the conflict.
3. Try to figure out if you can manually merge the two version --
maybe you have tried to fix the same thing, then decide which
solution is best, or maybe the conflict is really minor.
4. Fix it, remove the conflict markers and the text that is
not needed (remember to check if there are further conflicts
in the file)
5. Close the file, do `git commit -m "Fix merge conflict"` (without
a file name) and `git push`
6. If the conflicts are hard to resolve, maybe it is unclear what
you can remove or maybe very large parts of the file is affected,
don't hesitate to ask for help -- you might need to communicate
with the person who committed the text in conflict with yours or
you might want to consult someone more experienced.


## Create a new branch

If you want to introduce modifications or develop new features for a certain
application without having to change anything in the main branch `master`,
then you can create a new branch.


Let's say you are on the `master` branch and want to create a new branch.
Here's what you'll do:
```
git checkout -b <my branch name>.
```

This command will automatically create a new branch and then move you to
that branch. To confirm that your new branch was created you can use the
`git branch` command. The branch name with the asterisk next to it indicates
which branch you're pointed to at that given time.


Once you've created a branch, you can save your progress with `git add` and
`git commit` as usual. Remember that executing `git push` will upload the
changes _remotely_ to the current branch `<my branch name>` (_NB!_ You may
perhaps get a message from git and  be asked to push explicitly to the remote
version of the branch -- I usually just follow the syntax that is suggested
in that message).


Now if you are ready to make all the modifications in the branch to be part
of the main branch, then you can merge it to `master` by executing:

```
git checkout master
git merge <my branch name>
```

Note that at any point in time, you can switch to the `master` branch by
`git checkout master`. Similary, if you want to switch to a different branch
then do `git checkout <my branch name>`.

### Forking a repository

If you are working in a project with multiple users or in a project owned by
someone else, then it is useful to first fork the repository, that is, to make
your own copy of the repository that you can work with and do the changes you
need to do. In many projects, this is the only way of contributing changes
to a repository.

This is most easily done on repository web-page. For a GitHub repository,
you click on the `Fork` button at the top right of the repo page.
You will then be asked where to place it -- usually you want to place it in
your own workspace in an appropriate project. Then from your so created
*local* fork repository web page, you continue by cloning a _gwd_ and then
create a branch, implement, commit, and push changes, just as described above.
(If you want to share your changes with the repo owner, you can make a
_pull request_ -- currently not covered here.)

The original repository that you forked from is called the upstream
repository of your fork.

#### Updating a forked repository from the upstream repository

You might want to update your fork with any new changes in the upstream
repository. This is most easily done on the for repository web-page, by
clicking the `Fetch upstram` button at the upper right just above the
file list (Github); In Bitbucket, there will be a button
`Sync (X commits behind)` under repository details in the upper right
corner.

#### Removing a forked repository

Removing a forked repository will not affect the original repository at all.
Removal is best done on the webpage for the forked repository:

###### Github

1. Verify that you are in the right repository!
2. Click on the `Settings` tab near the top of the page.
3. Click `Options`in the left-hand menu.
4. Scroll down to the bottom of the pages.
5. In the `Fanger Zone`, click the `Delete this repository` button.
6. Follow the instructions to complete the deletion.

###### Bitbucket

1. Verify that you are in the right repository!
2. Click `Repository settings` in the left margin.
3. Click `Repository details` under `GENERAL` to the left of the Pages.
4. Click the rolol-down menu `manage repository` on the top right of the
page and click `Delete repository`.
5. Verify that you want to delete.


## Good `git` practices
It is good to keep your repo clean:

- Keep the code and the indata separated; **Never** put data on ’git’.
- Try to structure your files into logical subfolders, so that it is easy
to navigate the repo.
- Add a README file, preferably in _markdown_-format_ (.md_) in the repo
top level folder, that briefly explains the purpose and structure of the
repo. This README will display on the web-page of the repo.
- Also, add README files in strategic subfolders.

When developing code/text (i.e., changing files) in `git`, there
are some things one can do to avoid _merge conflicts_:
- Start each developing session by `git pull` so you don't work
on outdated files.
- End each session by `git push` for much the same reason.
- If you feel you need to do mayor changes that may take longer
time and maybe lots of try-outs, you should read up on
_branching in git_ and/or _forking repositories_ (above, or google it).

## Further reading
This crash course is intentionally rather NBIS lts project-specific.
There are several other, more general, crash courses or tutorials
that complements, improves or extends this text. You can find
them by googling.

Github provides a rather good and detailed `git`
[documentation](https://git-scm.com/docs), which can also be googled;
actually googling for the feature you are interested in is generally
a very productive idea.
Specifically, one can look at the commands `git log`, `git stash`,  
look more into _branching in git_, and the concept of _forking repositories_,
and how to make pull requests.

If you are setting up your own repo, to which others can contribute, you might
want to look up the _Continuous integration_ (_CI_) feature perovided by `git`
repo service providers. This can help you check for conflicts and run test on
incoming pull requests or commit and simplify maintaining your repository. But
this is definitely for a more advacned course:).


## Exercises

These will require that you have a Bitbucket user account.

1. Create a fork of the present repository (i.e., on the web page
https://github.com/NBISweden/knowledgetransfercrashcourses.git).  
Verify that the address of forked the repository and the original repository
are different.
2. Clone a git working directory (_gwd_) from your *fork* of the present
repository.  Can you see it?
3. `cd` into the _gwd_ and then into the _GitCrashCourse_ folder. Open the file
`README.md` in a text editor(e.g., atom or sublime). Change the title from  
_# Git Crash Course_  
to   
_# My Git Crash Course_  
Commit the changes to git with an appropriate message, and push them to the
forked repository. Update the *fork* webpage in your browser.  Can you see
the changes?
3. Create a branch called "_MyBranch_". Use a text-editor to create a file
`myfile` and add some text content. Add the file to `git`, commit and push it.
Update the fork page in the browser.  
Can you se the file?
6. Switch back to the `master` branch.
5. You will need the forked repository in the exercises of some of the other
crash courses. When you have completed them, you can delete the fork.
