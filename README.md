# README

## Setting up your working directory

To clone this repo run:
`git clone https://github.com/julia-neme/EAC_watermasses.git`

We have a .gitignore set up to ignore data files and figures, so in the EAC_watermass directory, create two additional directories called `results` and `data`. Your directory structure should look like this:
```
EAC_watermasses
  ├── data
  ├── docs
  ├── results
  └── src
```

## Create your branch (here you work)

`git checkout -b "your_name"`

And before working, always verify that you are on your own branch by running `git branch`. Once you are done working, add and commit the files (avoid using `git add --all` and add the files you want to push), push your branch to the repo:

`git push origin "your_name"`

On the github repository you should now see that a branch has pushed changes, and create a request to merge with the main branch.

## Sync your local branch to master

```
git checkout master
git pull origin master
git checkout "your_name"
git merge origin/master
```
