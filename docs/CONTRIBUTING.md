# Contributing to Biobase
Thank you for your interest in contributing to Biobase! We welcome contributions, whether you're fixing a bug, improving documentation, or adding new features. This guide will help you get started with contributing to the project.

## Table of Contents
- [Code of Conduct](#Code-of-Conduct)
- [How to Contribute](#How-to-Contribute)
- [Bug Reports](#Bug-Reports)
- [Feature Requests](#Feature-Requests)
- [Pull Requests](#Pull-Requests)
- [Code Style](#Code-Style)
- [License](#License)
- [Community](#Community)
## Code of Conduct
By participating in this project, you agree to abide by our [Code of Conduct](https://GitHub.com/lignum-vitae/biobase/blob/master/docs/CODE_OF_CONDUCT.md). Please take a moment to familiarize yourself with it.

## How to Contribute
### Bug Reports
If you find a bug or unexpected behavior, please open an issue in the GitHub repository. When reporting a bug, provide the following information:

- A description of the problem.
- Steps to reproduce the issue (if applicable).
- Any relevant error messages.
- The version of the library you're using.
- The Python version you're using.

### Feature Requests
If you have an idea for a new feature or enhancement, please open an issue describing the feature and why you think it would be useful. We encourage open discussions before starting to code a new feature.

### Pull Requests
To contribute code:

#### Note
Detailed below is the process of adding the repo as an upstream repo through your Command Line Interface (CLI).
However, GitHub allows you to sync your fork through their Web UI by navigating to the GitHub Page of your repo fork and clicking on the `Sync fork` button.
GitHub also has its own CLI that allows you to use the command `gh repo sync owner/cli-fork -b BRANCH-NAME`.
You can read more about that in the [GitHub Docs](https://docs.GitHub.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/syncing-a-fork) here.
If you go this route, you should be able to skip steps 4 and 9, as well as omit `upstream/main` from step 5 below. Push your changes directly to your GitHub fork.

#### 1. Open a new Issue following the above-mentioned guidelines.
#### 2. Fork the repository to your own GitHub account.
#### 3. Clone your fork locally:
```nginx
git clone https://GitHub.com/YOUR_USERNAME/biobase.git
```
#### 4. Keep your fork up to date with the main branch
```nginx
# Add the main Biobase repo as upstream
git remote add upstream git://GitHub.com/lignum-vitae/biobase.git
# Get latest changes
git fetch upstream
# Verify your remotes
git remote
```
#### 5. Create a new branch for your changes:
```nginx
git checkout -b feature-name upstream/main
```
#### 6. Make your changes in your local repository.
#### 7. Run your changes in your local environment.
```nginx
# Downloads project as editable, which allows local imports
python -m pip install -e .

# Example to run motif.py
python -m src.biobase.constants.analysis.motif
```

- Install the project as editable by going to the root directory
  and running the command `python -m pip install -e .` from the command line.
  Make sure you're in the folder that contains the pyproject.toml file.
- To run a specific file from the command line, run it using dot separators from the root of the project
  - For example, to run matrix.py, use the command `python -m src.biobase.matrix`
  - To run motif.py, use the command `python -m src.biobase.constants.analysis.motif`

#### 8. Commit your changes with a descriptive commit message:
```nginx
git add <filename>
git commit -m "Add feature: description of change"
```
#### 9. Rebase Your Development Branch on the Latest Upstream
```nginx
# Gets latest changes from main biobase project if you've set up an upstream branch as detailed above
git fetch upstream
# Make sure all is committed (or stashed) as necessary on this branch
git rebase -i upstream/main feature-name
```
You may need to resolve conflicts that occur when both a file on the development trunk and one of the files in your branch have been changed. 
Edit each conflicting file to resolve the differences, then continue the rebase. 
Each file will need to be "added" to mark the conflict as resolved:
```nginx
$ git add <filename>
$ git rebase --continue
```
#### 10. Push your branch to your fork on GitHub:
```nginx
git push -f origin feature-name
```
#### 11. Open a Pull Request (PR) from your branch to the main branch of the original Biobase repository on GitHub.
- You may need to click the `compare across forks` link under the `Compare changes` header that populates
  when you click `New pull request` to see your local repo fork.

#### 12. In your PR description, include:
- A summary of the changes.
- Any relevant issue numbers (e.g., fixes #123).
- Information about tests and validation.

## Code Style
We follow standard Python conventions (PEP 8) for code style. Some additional notes:

- Use meaningful variable and function names.
- Keep lines of code under 80 characters where possible.
- Make sure to update documentation if your changes affect the usage or API.
## License
By contributing to Biobase, you agree that your contributions will be licensed under the MIT License, as outlined in the [LICENSE](https://GitHub.com/lignum-vitae/biobase/blob/master/LICENSE) file.

## Community
We encourage contributions from everyone, and we strive to maintain a welcoming and inclusive community. If you have any questions, need help, or want to discuss ideas, feel free to reach out via issues or the repository discussions.

Thank you for contributing to Biobase! Your help improves the project for everyone!

