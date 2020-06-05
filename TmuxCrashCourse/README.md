# Tmux Crash Course

## Introduction
`tmux` is a virtual terminal program that can be used, for example, to keep 
long running jobs alive on UPPMAX and other clusters.

## Why?
A typical example is provided by running a `snakemake` workflow on th UPPMAX cluster _rackham_. When we start `snakemake`, a master process first starts. This master process handles all the admin of the workflow: it figures out which files needs to be produced,  which rules should be run, sends the specified runs as jobs to _rackham_'s cluster nodes, receives the results and evaluates what needs to be done next.

In other words, this master process needs to stay alive during the whole lifetime of the workflow, but most of the time all it does is waiting. So we don't want to send it as a job to a cluster node, but rather have it going on one of the _rackham_ login nodes.

However, processes started in the UPPMAX login terminal dies when we log out. We want to avoid that, but we usually can't stay logged in to UPPMAX from our laptop forever. Enter, a virtual terminal handler, like `tmux` or `screen`. We start a virtual terminal on the _rackham_ login node, start `snakemake` in the virtual terminal, and then _detach_ the virtual terminal. If we now log out from _rackham_, the virtual terminal stays alive and so does the `snakemake` process which will tick away while we do other things. We can log back in to the same login node on _rackham_ and _attach_ the virtual terminal and check the progress of `snakemake`.

## How?
There are a fairly few commands one need to know to run `tmux.`
Simply typing `tmux` will start a new virtual terminal. However, it is good practice to give each virtual terminal a name, in case we start more tmux-terminals and later need to identify which terminal is which.

### Using a better version of `tmux`
The `tmux` installed by default on uppmax is not among the latest versions, so it is good to get a more recent one using UPPMAX's module system, e.g.,
```
module load tmux/2.5
```

### Start a terminal
```
tmux new -s aGoodName
```
starts a new tmux terminal named _aGoodName_. Of course, you can change _aGoodName_ to any name you prefer.

### _Detach_ a terminal
To _detach_ the tmux terminal, hit
```
ctrl-b d
```
that is press the `ctrl` key and `b` at the same time, release and then press `d`. You're now back in your login terminal.

### Check what terminals you have created
```
tmux ls
```
This is good to recall names of terminals, among other things.

### _Attach_ a terminal
```
tmux a
```
_attaches_ the last used tmux terminal. If you want to _attach_ another terminal than the last one, you need to provide its name
```
tmux a -t aGoodName
```

### Stop/kill a terminal
If you have attached the terminal, you just type
```
exit
```
and the tmux terminal cease to be. Sometimes you want to kill a terminal while on the login terminal:
```
tmux kill-session aGoodName
```
fixes that.

### Scrolling
An annoying thing with `tmux`is that you can't scroll up to look at lines higher in the terminal window. However, there's a solution for that (if you use the newer version available from UPPMAX module system, at least). Typing
```
ctr-b page-up
```
(i.e., press the `ctrl` key and `b` followed by the `page-up` key or the corresponding key combination -- on my mac it is `shift` and the `up-arrow`) allows you to scroll up. You must then tell `tmux` that you quit scrolling by simply typing
```
q
```

I think that's most of the things you ned to know. Good luck
