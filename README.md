# Data Repository for the Code Comparison


This is the repository of a collaboration of several supernova radiative transfer codes. 
This repository will be updated from time to time. There is a code called plot_XXX that can be used to 
easily make comparison plots with the following instructions.

```bash

./plot_XXX toy01 lbol

```

[![Build Status](https://dev.azure.com/sn-rad-trans/data1/_apis/build/status/sn-rad-trans.data1?branchName=master)](https://dev.azure.com/sn-rad-trans/data1/_build/latest?definitionId=2&branchName=master)

# Updating the Data

Follow the steps below if you would like to update the data yourself using github.

1. Fork the Data1 repository (using the fork button on the github website)  -- This creates a copy of the Data1 repository on your github account.
2. `cd` into your forked Data1 repo.
3. `git remote add upstream git://github.com/sn-rad-trans/data1` -- This connects your forked version of the Data1 repo to the original version.
4. `git pull upstream master` -- This makes sure that your forked version is up to date with the original repo.
5. `git branch data_update` -- This creates a new branch of your forked repo named "data_update." Consider adding the name of your code to the name to make the branch name as descriptive as possible.
6. `git checkout data_update` -- This checks out your new branch, so changes you make will only affect the new branch of your forked repo. (If you break things beyond repair, just delete the branch, checkout your master branch, and start again.)
7. Make your changes. Replace the old files with the updated version of the files.
8. `git status` -- Shows you what you have changed in the repo.  
9. `git add <all new files>` -- You need to run `git add` with each new file listed under the "Untracked Files" category of the `git status` report.
10. `git commit -am "<Short, descriptive update message.>"` -- This commit all changes you have made.  Please replace <Short, descriptive update message> with something useful like: "Updated Stella files."
11. `git push --set-upstream origin data_update_stella` -- This will push your changed branch to your forked repo.
12. Navigate to your forked version of the repo. Under the "Branch" dropdown menu select your new branch "data_update"
13. Click "New Pull Request" button, add a brief description of the changes you have made, and click "Create Pull Request" -- This will notify the team managing the original repo that you have made changes that are ready to merge.
