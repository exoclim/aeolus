Contributor's Guide
===================

Aeolus takes contributions using pull requests on `GitHub <https://github.com/exoclim/aeolus/pulls>`_.

Setting up repositories
-----------------------
0. Create an account on GitHub.
1. Fork aeolus on GitHub, so a copy of it appears in your GitHub profile and has a remote URL like `https://github.com/your_user_name/aeolus <https://github.com/your_user_name/aeolus>`_.
2. Git-clone it to your computer using the new URL (or `change the remote <https://docs.github.com/en/github/getting-started-with-github/getting-started-with-git/managing-remote-repositories#changing-a-remote-repositorys-url>`_ of your existing folder).


Developing the code
-------------------
0. It is often more convenient to install the library you are developing in the "developer" mode. To do this, navigate to the cloned copy of the repository and type :code:`pip install -e .` (Make sure you are using :code:`pip` in the correct python environment and note the dot at the end!) Alternatively, run :code:`python setup.py develop` (no dot at the end).
1. It is also better to create a new branch instead of working on the main branch (e.g. :code:`git checkout -b new_branch_for_cool_addition`). If your contribution is relatively small, such as fixing one small bug or typo, you can skip this step.
2. Make changes, commit them and push to *your* remote. (Which again, should have your URL. You can check it by typing :code:`git remote -v`).
3. Test the changes (see below)!
4. Go to your aeolus page on GitHub and create a pull request (a prompt for this should show up on top of the page).
5. Wait for the automatic tests to pass and for aeolus developers to review your pull request. If any further changes are needed, make the necessary changes and commit to your forked remote repository. The pull request will update automatically.


Checking the code formatting
----------------------------
0. Install :code:`pre-commit` in the current environment: :code:`pip install pre-commit`.
1. Check the formatting style by running :code:`pre-commit` on all files: :code:`pre-commit run --all-files`.
2. The errors should be corrected automatically, otherwise fix them manually.
3. On your next commit, the :code:`pre-commit` hook should work automatically.


Testing the code
----------------
0. Run your own tests to make sure the code does what you expect (make sure you are importing the "right" version of aeolus).
1. Download the test data (from the top level of the aeolus directory) :code:`git clone https://github.com/exoclim/aeolus_data.git tests/data`
2. Run aeolus tests :code:`pytest -c setup.cfg --cov=aeolus --cov-report term-missing -v`
3. Fix all the errors, if there are any!.


Contributing to the documentation
---------------------------------
We welcome contributions to improve the documentation of aeolus.

0. Please make yourself familiar with the `reST formatting guide <https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html>`_.
1. Follow the steps above to create a fork on GitHub.
2. Make changes. Do not forget to add any relevant references to papers or textbooks.
3. Build the documentation locally.
    a. Create a separate conda environment :code:`conda env create -f docs/environment.yml`
    b. :code:`cd` to the `docs/` directory.
    c. Build the documentation :code:`make html`
    d. :code:`cd` to the `_build/html` sub-directory.
    e. Start an HTTP server using a random port :code:`python -m http.server 1234`
    f. Open the :code:`http://localhost:1234` page in the browser. Use CTRL-F5 to refresh the page if needed.
    g. To close the session, use CTRL-C in the command line where the HTTP server was started.
4. Create a pull request and wait for the aeolus team to review it.
