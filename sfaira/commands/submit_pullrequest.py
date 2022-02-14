import logging
import subprocess
import os
import shutil
import sys

from sfaira.commands.questionary import sfaira_questionary
from rich import print


log = logging.getLogger(__name__)


class PullRequestHandler:

    def __init__(self, path_loader):
        self.WD = os.path.dirname(__file__)
        self.path_loader = path_loader
        self.loader_name_list = []
        self.loader_name = None
        self.gh_token = ""

    def submit_pr(self):
        self._check_container_and_loaders()
        self._gh_authenticate()
        self._clone_sfaira_and_move_dataloader()
        self._submit_pr()

    def _check_container_and_loaders(self):
        sfaira_container = os.getenv('SFAIRA_DOCKER')
        if sfaira_container and self.path_loader and os.path.isdir(self.path_loader):
            for d in os.listdir(self.path_loader):
                if os.path.isdir(os.path.join(self.path_loader, d)) and d.startswith(("d10_", "dno_doi_")):
                    self.loader_name_list.append(d)
            if not self.loader_name_list:
                print("[bold red]No data loaders found in loader directory. Aborting.")
                sys.exit()
        else:
            print('[bold red]This function can only be called when running inside the sfaira CLI docker container. '
                  'Aborting.')
            sys.exit()

    def _gh_authenticate(self) -> None:
        """
        Guides the user to authenticate with the GitHub CLI
        """
        # Skip login procedure if already logged in
        if subprocess.run(["gh", "auth", "status"], check=False, text=True, shell=False).returncode == 0:
            subprocess.run(["gh", "auth", "setup-git"], check=True, text=True, shell=False)
            print("[bold green]Already authenticated with GitHub.")
            return
        print(
            "[bold blue]Please authenticate with GitHub. Hint: you can generate a Personal Access Token here: "
            "https://github.com/settings/tokens\nThe minimum required scopes are 'repo', 'read:org', 'workflow'."
        )
        returncode = 1
        while returncode != 0:
            gh_token = sfaira_questionary(function='password',
                                          question='Please enter your GitHub token. '
                                                   'Leave blank to interactively authenticate with GitHub.',
                                          default='')
            if gh_token == "":
                gh_call = subprocess.run(
                    ["gh", "auth", "login", "--hostname", "github.com"], check=False, text=True, shell=False)
            else:
                gh_call = subprocess.run(
                    f"echo {gh_token} | gh auth login --with-token", check=False, text=True, shell=True)
            returncode = gh_call.returncode
        self.gh_token = gh_token
        subprocess.run(["gh", "auth", "setup-git"], check=True, text=True, shell=False)
        print("[bold green]Successfully authenticated with GitHub.")

    def _clone_sfaira_and_move_dataloader(self) -> None:
        """
        Clones the sfaira repo, creates a new branch and moves dataloader into the right location.
        """
        # Clone sfaira
        subprocess.run(["rm", "-rf", "/root/sfaira"], check=False, text=True, shell=False)
        subprocess.run(
            ["gh", "repo", "clone", "theislab/sfaira", "/root/sfaira/"], check=True, text=True, shell=False)
        # Get loader name
        if len(self.loader_name_list) == 1:
            self.loader_name = self.loader_name_list[0]
        else:
            self.loader_name = sfaira_questionary(
                function='select',
                question='Multiple dataloaders detected in the loader directory. '
                         'Which one do you want to submit as a pull request?',
                choices=self.loader_name_list
            )
        # Create new branch in sfaira git repo
        subprocess.run(
            ["git", "checkout", "-b", f"dataset/{self.loader_name}"],
            check=True, text=True, shell=False, cwd="/root/sfaira/"
        )
        # Copy loader
        shutil.copytree(
            src=os.path.join(self.path_loader, self.loader_name),
            dst=os.path.join("/root/sfaira/data/dataloaders/loaders", self.loader_name)
        )

    def _submit_pr(self):
        """
        Adds and commits the dataloader in git and creates a pull request
        """
        # Make sure git credentials are set and get them from github api if not
        current_email = subprocess.run(["git", "config", "user.email"], check=False, text=True, shell=False,
                                       stdout=subprocess.PIPE).stdout
        if not current_email:
            git_email = subprocess.run(["gh", "api", "/user/public_emails", "-q", ".[0].email"], check=True, text=True,
                                       shell=False, stdout=subprocess.PIPE).stdout.strip()
            subprocess.run(["git", "config", "--global", "user.email", git_email], check=True, text=True, shell=False)
        current_user = subprocess.run(["git", "config", "user.name"], check=False, text=True, shell=False,
                                      stdout=subprocess.PIPE).stdout
        if not current_user:
            git_user = subprocess.run(["gh", "api", "/user", "-q", ".login"], check=True, text=True, shell=False,
                                      stdout=subprocess.PIPE).stdout.strip()
            subprocess.run(["git", "config", "--global", "user.name", git_user], check=True, text=True, shell=False)
        # Add and commit the dataloader
        subprocess.run(["git", "add", "*"], check=True, text=True, shell=False, cwd="/root/sfaira/")
        subprocess.run(["git", "commit", "-m", f"[from sfaira cli] add dataloader {self.loader_name}"],
                       check=True, text=True, shell=False, cwd="/root/sfaira/")
        # Create pullrequest (authenticate again beforehand if gh_token was provided before)
        create_pr_str = f"gh pr create --base dev --title {self.loader_name} " \
                        f"--body 'This PR was created by the sfaira CLI adding dataset {self.loader_name}' " \
                        f"--label dataset"
        if self.gh_token != "":
            subprocess.run(f"echo {self.gh_token} | gh auth login --with-token && {create_pr_str}",
                           check=True, text=True, shell=True, cwd="/root/sfaira/")
        else:
            subprocess.run(create_pr_str, check=True, text=True, shell=True, cwd="/root/sfaira/")
        subprocess.run(
            ["rm", "-rf", os.path.join(self.path_loader, self.loader_name)], check=False, text=True, shell=False)
        print("[bold green]Your PR was successfully submitted. Feel free to add further comments to it using the URL "
              "in the line above.")
