# Working with Git
### What is Git
- A version control system that saves a snapshot of your project every time you make changes
- Each change has to be documented with a description of the changes
- If multiple people are working on the same project (or if one person works on it in multiple places), changes can be synced
- Multiple versions of can be maintained and changes can be merged across versions
- Syncing requires storing the project on a Git server (eg, [GitHub](https://github.com/) or [GitLab](https://about.gitlab.com/))
- Since code is stored online, it's easy to share with others
- GitHub is a free service, but normally requires that all code to be publicly available (private repos are available for academic use)
    - this isn't normally a problem, though you'll want to ensure that your **data** are not included in the repository

### [The preferred workflow for Git](http://nvie.com/posts/a-successful-git-branching-model)

![GitFlow Branching Model](http://nvie.com/img/git-model@2x.png)

### [A realistic workflow for Git](https://xkcd.com/1597)

![XKCD Branching Model](http://imgs.xkcd.com/comics/git.png)

### Working with GitHub Desktop
- Installing GitHub Desktop:
    - Available from https://desktop.github.com
    - Create [GitHub](github.com) account and log into GitHub Desktop using these credentials
- Adding an Existing Project:
    - Click on project selector in upper right corner of GitHub Desktop and select ```New Project```
    - Click on ```Existing Directory``` and navigate to folder, then click ```Select Folder``` and ```Create Project```
    - Click on ```+``` in upper left corner of GitHub Desktop, then ```Create```
    - Navigate to folder and click ```OK```
    - Give repository a name and create (note that redundant names in different locations on your file system will clash on the server)
- Starting from Scratch:
    - Click on project selector in upper right corner of GitHub Desktop and select ```New Project```
    - Click on ```New Directory```, then ```Empty Project```, then pick a project name and location
    - Click on ```+``` in upper left corner of GitHub Desktop, then ```Create```
    - Navigate to folder and click ```OK```
    - Give repository a name and create
- To make changes to the repo:
    - Click on ```Changes``` at top of window
    - Select files to add
    - Write summary of changes and (optional) description
    - Click ```commit to master```
- To revert changes:
    - Click on ```History``` at top of window
    - Navigate through commits to find the offending edit (deletions are marked red, additions are green)
    - Click on ```Revert```
- To sync local copy to remote server:
    - Click ```Sync```

### Resolving conflicts
- If sync fails because changes on remote server can't be automatically merged with local changes, open file in your text editor and edit conflicts manually
- Conflicts look like this:

    ```
    <<<<<<<<<<<< HEAD
    local version of code
    =============
    remote version
    >>>>>>>>>>>> origin/master
    ```

- Edit file to include whatever parts from the local and remote versions you please (and remove markers added by Git)
- Return to GitHub Desktop, commit changes, and redo sync
- If all else fails, save local files in another location, delete the project, clone fresh copy from server, then add back in your local changes
    
### [Working with branches](https://git-scm.com/book/en/v2/Git-Branching-Basic-Branching-and-Merging)
- I usually don't bother creating additional branches beyond the local and remote master versions
- If multiple people are working on the same project, it's nice if they can work on the same files independently without having to merge everything together each time they make a commit
    - This can acomplished by creating seperate branches for each person, then merging everything together once everyone has completed their part of the project
- Can also be useful to "snapshot" your project, eg; when submitting for publication
- Create a new branch by clicking on the ```Create new branch``` is in the upper left, next to the ```master``` dropdown menu
- Give your branch a name and click ```OK```
- Any changes will now be confined to your new branch
- The ```Sync``` button will now be replaced by a ```Publish``` button because there is no copy of your new branch on the server to sync to
- After clicking on ```Publish```, the ```Sync``` button will come back and any changes made to your files will not conflict with changes on the master version on the server
- You can now use the drop down menu in the upper left to switch between branches, which will automatically change the code in your folder
- By some black magic, this changes the files that are available in the finder
- When you are ready to merge branches back together, return to master, click on ```compare``` then select the branch you want to merge from, and click on ```Update from ...```
- It's usually a good idea to merge from master to your new branch and resolve any conflicts, then merge back into master

### Using git from the command line 
- Copy project url and clone locally:

    ```git clone https://github.com/hobrien/RNAseqTools.git```

- List modified files:

    ```git status -s```

- Add changes to project:

    ```git add FILENAME```
    ```git commit -m "initial commit"```

    - can also commit without -m flag, which opens the commit message in your text editor of choice
- Update local version with remote changes:

    ```git pull```

- Update remote repo with local committed changes:

    ```git push```
    - always sync in this order

- Undo changes to a file:

    ```git checkout -- FILENAME```

- Create new branch:
    - ```git branch -b Version2```
- Switch branches:
    - ```git checkout master```
- Merge branches:
    - ```git merge Version2```
    - ```git branch -d Version2```

### A few additional considerations
- Bookmarking snapshots (example: the version used for the analyses in a paper submission)
    - Make zipped copy of directory
    - Make a note the number of the last commit before submission (eg; ce289b3)
    - Create a separate branch for each submission
    - Tag commit
        - Go online and click on the "0 releases" button beneath the project description
        - Click on "Create a new release"
        - Give release a name (unfortunately only names like "v1.0" are supported)
        - Click on "This is a pre-release" if it's not the final version of the paper, then on "Publish release"
        - This will automatically create a ziped version of all the code that can be navigated to by clicking on the "1 release" button for the project
- Excluding private files
    - When you initialise a repo, a .gitignore file will automatically be created with a list of files that will never be synced to the server
    - Files on the list will not show up on the lists of modified files to add to a commit
    - It's a good idea to list any data files that you don't want to make public in this file, along with any config files that contain things like login credentials
    - It is a *really* good idea to keep this list up to date, because it's easy to accidentally commit all untracked files in your project folder
        - If you do this with genome-scale files, it will break *everything*
        - It's also not trivial to delete files from [github](github.com) once you've committed them, because git wants to keep a record of all changes made to all files
    - Of course, any files that are not part of the repo will need to be synced some other way
- README, LICENSE
    - When a repo is created online, you have the option of creating a README.md file that is displayed automatically on the landing page for your project
    - If you are creating a repo locally, you'd need to do this manually with a text editor
    - You are also given the option of adding a license to the repo, though the available options are tailored to software (I would suggest using the Creative Commons [CC-BY](https://creativecommons.org/licenses/by/2.0/uk/legalcode) license for most biology projects)
- Meta-notebooks
    - Notes can be organised in a logical order, without interleaving unrelated work that was done on the same day
    - Commit history preserves a chronological record
    - Additional explanation, comments, etc can be recorded in commit message without cluttering notes
    - If you're fortunate enough to have lab mentors, they can easily review your notes and leave comments
    - If using the free version of GitHub, you may need to be a bit coy about exactly what you are working on to prevent your competetors from stealing you work
- Using Git for manuscript writing
    - Can be used to avoid this:
    
    ![tweet](https://raw.githubusercontent.com/MixedModels/LearningMLwinN/master/ScreenShots/tweet.png)
    - Merging changes by different co-authors is done automatically
    - Keeps a record of all changes, who made them, and why
    - Discussions about particular changes can be had in the comment thread for that change instead of in email threads
    - Issue tracker can be used as a to-do list, and also includes comment threads for discussion of proposed tasks
    - Branching can be used to allow a co-author to make a series of commits that can all be merged at once
    - This may well be something you'd want a private repo for

<br>
<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br /><span xmlns:dct="http://purl.org/dc/terms/" href="http://purl.org/dc/dcmitype/Text" property="dct:title" rel="dct:type">Working with Git</span> by <span xmlns:cc="http://creativecommons.org/ns#" property="cc:attributionName">Heath O'Brien</span> is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.