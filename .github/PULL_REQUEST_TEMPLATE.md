## Overview

Short description of what the PR addresses, eg. new feature or bug fix.

## PR checklist
 - [ ] Adds documentation for the new scripts, functions you added
 - [ ] If scripts are inccluded in the docker container, make sure to build and push to DockerHub
 - [ ] Link issues that this PR addresses from the `Linked issues` section in the bottom right of the Pull Request page (after the PR is created)
 - [ ] Rebuild and push the container if you have updated the Dockerfile

## My changes are final, who do I tell?

Once your changes are final, request a PR review from the right side of the PR page, from [@cgpu](https://github.com/cgpu) or another relevant person from the team depending on the pipeline.

----
## My PR is approved, now what?

If your PR has been approved, here are the last steps to finalize your contribution
- Do not add any commits after the approved review. 
- Need to modify further after your PR has been approved? Make sure to add a note so that the reviewer and other collaborators can see.
- Merge your changes: at the bottom of the PR page, switch to `Squash and Merge`. 
**NOTE**: Avoid selecting `Create a merge commit`. This introduces an extra commit with an uninformative automatic message: _"Merge branch **this** of https://github.com/lifebit-ai/<repo name> into **that**"_. We want to avoid this to make sure the commit history stays informative.