# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SPHINXPROJ    = MultiphaseFlowSimulator
SOURCEDIR     = .
BUILDDIR      = _build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
html: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

# default value: deploy_gh_pages
ifndef DEPLOY_DEFAULT
DEPLOY_DEFAULT = deploy_gh_pages
endif

# The deployment directory to be deployed
ifndef DEPLOY_DIR
DEPLOY_DIR      = _deploy
endif

# The heroku deployment directory to be deployed
# we must create this separated dir to avoid any conflict with _deploy (rsync and gh_pages)
ifndef DEPLOY_DIR_HEROKU
DEPLOY_DIR_HEROKU = _deploy_heroku
endif

# Copy contents from $(BUILDDIR) to $(DEPLOY_DIR)/$(DEPLOY_HTML_DIR) directory
ifndef DEPLOY_HTML_DIR
DEPLOY_HTML_DIR = docs
endif


## -- Rsync Deploy config -- ##
# Be sure your public key is listed in your server's ~/.ssh/authorized_keys file
ifndef SSH_USER
SSH_USER       = user@domain.com
endif

ifndef SSH_PORT
SSH_PORT       = 22
endif

ifndef DOCUMENT_ROOT
DOCUMENT_ROOT  = ~/website.com/
endif

#If you choose to delete on sync, rsync will create a 1:1 match
ifndef RSYNC_DELETE
RSYNC_DELETE   = false
endif

# Any extra arguments to pass to rsync
ifndef RSYNC_ARGS
RSYNC_ARGS     =
endif

## -- Github Pages Deploy config -- ##

# Configure the right deployment branch
ifndef DEPLOY_BRANCH_GITHUB
DEPLOY_BRANCH_GITHUB = gh-pages
endif

#if REPO_URL_GITHUB was NOT defined by travis-ci
ifndef REPO_URL_GITHUB
# Configure your right github project repo
REPO_URL_GITHUB       = git@github.com:victoraraujo01/multiphase-flow-simulator.git
endif

## -- Heroku Deployment Config -- ##

ifndef REPO_URL_HEROKU
# Configure your right heroku repo
# REPO_URL_HEROKU = git@heroku.com:spxd.git
endif


## end deployment configuration, don't edit anything below this line ##
#######################################################################

ifeq ($(RSYNC_DELETE), true)
RSYNC_DELETE_OPT = --delete
endif

init_gh_pages:
	@rm -rf $(DEPLOY_DIR)
	@mkdir -p $(DEPLOY_DIR)
	@cd $(DEPLOY_DIR); git init;\
		echo 'sphinx docs comming soon...' > index.html;\
		touch .nojekyll;\
		git add .; git commit -m "sphinx docs init";\
		git branch -m $(DEPLOY_BRANCH_GITHUB);\
		git remote add origin $(REPO_URL_GITHUB);
	@cd $(DEPLOY_DIR);\
		if ! git ls-remote origin $(DEPLOY_BRANCH_GITHUB) | grep $(DEPLOY_BRANCH_GITHUB) ; then \
			echo "Preparing Github deployment branch: $(DEPLOY_BRANCH_GITHUB) for the first time only...";\
			git push -u origin $(DEPLOY_BRANCH_GITHUB);\
		fi

setup_gh_pages: init_gh_pages
	@echo "Setting up gh-pages deployment..."
	@cd $(DEPLOY_DIR);\
		git fetch origin;\
		git reset --hard origin/$(DEPLOY_BRANCH_GITHUB);\
		git branch --set-upstream $(DEPLOY_BRANCH_GITHUB) origin/$(DEPLOY_BRANCH_GITHUB)
	@echo "Now you can deploy to Github Pages with 'make generate' and then 'make deploy_gh_pages'"

init_heroku:
	@rm -rf $(DEPLOY_DIR_HEROKU)
	@mkdir -p $(DEPLOY_DIR_HEROKU)
	@cd $(DEPLOY_DIR_HEROKU); git init;\
		cp -r ../.deploy_heroku/* .;\
		echo 'sphinx docs comming soon...' > public/index.html;\
		git add .; git commit -m "sphinx docs init";\
		git remote add origin $(REPO_URL_HEROKU);
	@cd $(DEPLOY_DIR_HEROKU);\
		if ! git ls-remote origin master | grep master ; then\
			echo "Preparing Heroku deployment for the first time only...";\
			git push -u origin master;\
		fi

setup_heroku: init_heroku
	@echo "setting up heroku deployment..."
	@cd $(DEPLOY_DIR_HEROKU);\
		git fetch origin;\
		git reset --hard origin/master;\
		git branch --set-upstream master origin/master
	@echo "Now you can deploy to Heroku with 'make generate' and then 'make deploy_heroku'"

generate: html

prepare_rsync_deployment:
	@echo "Preparing rsync deployment..."
	@mkdir -p $(DEPLOY_DIR)/$(DEPLOY_HTML_DIR)
	@echo "Copying files from '$(BUILDDIR)/html/.' to '$(DEPLOY_DIR)/$(DEPLOY_HTML_DIR)'"
	@cp -r $(BUILDDIR)/html/. $(DEPLOY_DIR)/$(DEPLOY_HTML_DIR)

deploy_rsync: prepare_rsync_deployment
	@echo "Deploying on rsync now..."
	rsync -avze 'ssh -p $(SSH_PORT)' --exclude-from $(realpath ./rsync_exclude) $(RSYNC_ARGS) $(RSYNC_DELETE_OPT) ${DEPLOY_DIR}/ $(SSH_USER):$(DOCUMENT_ROOT)

prepare_gh_pages_deployment:
	@echo "Preparing gh_pages deployment..."
	@echo "Pulling any updates from Github Pages..."
	@cd $(DEPLOY_DIR); git pull;
	@mkdir -p $(DEPLOY_DIR)/$(DEPLOY_HTML_DIR)
	@echo "Copying files from '$(BUILDDIR)/html/.' to '$(DEPLOY_DIR)/$(DEPLOY_HTML_DIR)'"
	@cp -r $(BUILDDIR)/html/. $(DEPLOY_DIR)/$(DEPLOY_HTML_DIR)

deploy_gh_pages: prepare_gh_pages_deployment
	@echo "Deploying on github pages now..."
	@cd $(DEPLOY_DIR); git add -A; git commit -m "docs updated at `date -u`";\
		git push origin $(DEPLOY_BRANCH) --quiet
	@echo "Github Pages deploy was completed at `date -u`"

prepare_heroku_deployment:
	@echo "Preparing heroku deployment..."
	@echo "Pulling any updates from Heroku..."
	@cd $(DEPLOY_DIR_HEROKU); git pull;
	@mkdir -p $(DEPLOY_DIR_HEROKU)/public/$(DEPLOY_HTML_DIR)
	@echo "Copying files from .deploy_heroku to $(DEPLOY_DIR_HEROKU)"
	@cp -r .deploy_heroku/. $(DEPLOY_DIR_HEROKU)
	@echo "Copying files from '$(BUILDDIR)/html/.' to '$(DEPLOY_DIR_HEROKU)/public/$(DEPLOY_HTML_DIR)'"
	@cp -r $(BUILDDIR)/html/. $(DEPLOY_DIR_HEROKU)/public/$(DEPLOY_HTML_DIR)


deploy_heroku: prepare_heroku_deployment
	@echo "Deploying on heroku now..."
	@cd $(DEPLOY_DIR_HEROKU); git add -A; git commit -m "docs updated at `date -u`";\
		git push origin master --quiet
	@echo "Heroku deployment was completed at `date -u`"


deploy: $(DEPLOY_DEFAULT)

gen_deploy: generate deploy
