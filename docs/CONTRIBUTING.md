# Contributing to Biobase

Thank you for your interest in contributing to Biobase! We welcome contributions,
whether you're fixing a bug, improving documentation, or adding new features.
This guide will help you get started with contributing to the project.

We recommend using [uv](https://docs.astral.sh/uv/getting-started/installation/)
for managing dependencies and running the project locally.

After installing uv, check out [uv commands](https://docs.astral.sh/uv/reference/cli/).
The most important will be `uv run` to run project files, `uv sync` to sync dependencies,
`uv add` to add dependencies, `uv remove` to remove dependencies, and `uv pip`.

## Table of Contents

- [Code of Conduct](#code-of-conduct)
- [How to Contribute](#how-to-contribute)
- [Bug Reports](#bug-reports)
- [Feature Requests](#feature-requests)
- [Pull Requests](#pull-requests)
- [Code Style](#code-style)
- [License](#license)
- [Community](#community)

## Code of Conduct

By participating in this project, you agree to abide by our
[Code of Conduct](https://github.com/lignum-vitae/biobase/blob/main/docs/CODE_OF_CONDUCT.md).
Please take a moment to familiarize yourself with it.

## How to Contribute

### Bug Reports

If you find a bug or unexpected behavior, please open an issue in the GitHub
repository. When reporting a bug, provide the following information:

- A description of the problem.
- Steps to reproduce the issue (if applicable).
- Any relevant error messages.
- The version of the library you're using.
- The Python version you're using.

### Feature Requests

If you have an idea for a new feature or enhancement, please open an issue
describing the feature and why you think it would be useful.
We encourage open discussions before starting to code a new feature.

### Pull Requests

To contribute code:

#### Note

Detailed below is the process of adding the repo as an upstream repo through your
Command Line Interface (CLI).
However, GitHub allows you to sync your fork through their Web UI by navigating
to the GitHub Page of your repo fork and clicking on the `Sync fork` button.
GitHub also has its own CLI that allows you to use the command
`gh repo sync owner/cli-fork -b BRANCH-NAME`.
You can read more about that in the
[GitHub Docs](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/syncing-a-fork)
here.
If you go this route, you should be able to skip steps 4 and 9, as well as omit
`upstream/main` from step 5 below. Push your changes directly to your GitHub fork.

#### 1. Open a new Issue following the above-mentioned guidelines

#### 2. Fork the repository to your own GitHub account

#### 3. Clone your fork locally:

```nginx
git clone https://GitHub.com/YOUR_USERNAME/biobase.git
```

#### 4. Keep your fork up to date with the main branch

```nginx
# Add the main Biobase repo as upstream
git remote add upstream https://github.com/lignum-vitae/biobase.git
# Get latest changes
git fetch upstream
# Verify your remotes
git remote
```

#### 5. Create a new branch for your changes:

##### Choose ONE of the following commands

```nginx
# Creates a new branch that stays in sync with the main repository
git checkout -b feature-name upstream/main

# Checks out existing branch if you already have a branch locally
git checkout feature-name
```

#### 6. Make your changes in your local repository

#### 7. Run your changes in your local environment

```nginx
# Downloads project as editable, which allows local imports
# This command must be run from the root/same level as the pyproject.toml
# [dev] also installs the optional dev dependencies
uv pip install --editable ".[dev]"

# Ensure dependencies match the project's using uv sync
uv sync --all-extras

# Example to run motif.py
pwd # src
cd biobase/analysis
uv run motif.py
```

- Install the project as editable by going to the root directory
  and running the command `uv pip install --editable .` from the command line.
  Make sure you're in the folder that contains the pyproject.toml file.
- To run a specific file from the command line, run it using uv and the appropriate file path
  - For example, to run sub_matrix.py if you are currently in the matrix folder,
  use the command `uv run sub_matrix.py`
  - To run motif.py, use the command `uv run motif.py`

#### 8. Commit your changes with a descriptive commit message:

```nginx
# Gets latest changes from main biobase project if you've set up an upstream branch as detailed above
git fetch upstream
# We recommend individually adding each file with modifications
git add <filename>
# Commit files after all files with modifications have been added
git commit -m "Add feature: description of change"
```

#### 🚨 Using "git add ." when staging changes

While `git add .` is convenient for adding all modified files, it can lead to
messy commits. Consider using it only when:

- You've reviewed all changes
- You're certain about each modification
- You've checked git status first
- Your .gitignore is properly configured

#### 9. Rebase Your Development Branch on the Latest Upstream

```nginx
# Make sure all is committed (or stashed) as necessary on this branch
git rebase -i upstream/main feature-name
```

You may need to resolve conflicts that occur when both a file on the development
trunk and one of the files in your branch have been changed.
Edit each conflicting file to resolve the differences, then continue the rebase.
Each file will need to be "added" to mark the conflict as resolved:

```nginx
# Resolve conflicts in each file, then:
git add <resolved-filename>
git rebase --continue
```

#### 10. Push your branch to your fork on GitHub:

```nginx
git push -f origin feature-name
```

#### 11. Open a Pull Request (PR) from your branch to the main branch of the original Biobase repository on GitHub

- You may need to click the `compare across forks` link under the `Compare changes`
header that populates when you click `New pull request` to see your local repo fork.
- Remember to add tests to the tests directory for any new functionality and ensure
that all tests pass when you run `pytest tests` from the root directory

#### 12. In your PR description, include:

- A summary of the changes.
- Any relevant issue numbers (e.g., fixes #123).
- Information about tests and validation.

## Code Style

We follow [Black](https://black.readthedocs.io/en/stable/getting_started.html)
for code formatting.
Black enforces a consistent style automatically, so please run it before committing.

Some additional notes:

- Use meaningful variable and function names.
- Keep lines readable (Black defaults to 88 characters).
- Make sure to update documentation if your changes affect the usage or API.

You can format the codebase by running this from the root:

```bash
black .
```

## License

By contributing to Biobase, you agree that your contributions will be licensed
under the MIT License, as outlined in the [LICENSE](https://github.com/lignum-vitae/biobase/blob/main/LICENSE)
file.

## Community

We encourage contributions from everyone, and we strive to maintain a welcoming
and inclusive community. If you have any questions, need help, or want to discuss
ideas, feel free to reach out via issues or the repository discussions.

Thank you for contributing to Biobase! Your help improves the project for everyone!
