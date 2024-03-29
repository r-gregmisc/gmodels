# gmodels Version 2.19.0  - 2024-04-04

New features:
  - Expose (and fix) `est_p_ci` which formats the estimate, p-value, and confidence interval into text.

Other Changes:
  - Upgrade to current package development tooling, including:
  - Use `roxygen2` for documentation and NAMESPACE
  - Use `NEWS.md` instead of `NEWS`
  - Add additional contributors to authors list and use new format.

# gmodels Version 2.18.1 - 2018-06-25

Other Changes:
  - Remove soft-links for NEWS and ChangeLog files for platform portability.

# gmodels Version 2.18.0 - 2018-06-19

Bug fixes:
  
  - ci.binom() was using an incorrect method for calcuating binomial
    conficence intervals.  It now calculates the Clopper-Pearson 'exect'
    interval, which is *conservative* due to the discrete nature of the
    binomial distribution.

Other Changes:
  
  - Support for lme4 objects has been removed due to incompatible
    changes to the lme4 package.

# gmodels Version 2.16.0 - 2014-07-24

New features:
  
  - The estimable() function now returns objects that are of class
    'estimable'.
  
  - The confidence interval function ci() now has a method for
    'estimable' objects, with the same layout as for 'lm' objects,
    making it easier to combine confidence information about model
    parameters and estimable functions into a single table.


# gmodels Version 2.15.5 - 2013-07-18

Bug fixes:

- Correct error in estimable.mlm() that caused it to always fail.  Added
  test code to prevent future issues.

Other Changes:

- Update man page file for ci() to current Rd syntax.
- Remove unused argument to ci.mer()


# gmodels Version 2.15.3 - 2012-06-27

Bug fixes:

- Update est.mer() to work with "mer" object changes introduced in
  lme4 version 0.999999-0.


# gmodels Version 2.15.2 - 2012-04-19

Bug fixes:

- Update est.mer() to work with recent versions of lme4  which changed
  'mer' objects from S3 to S4 class

- Changes to pass new R CMD check tests

- The 'Design' package has been replaced my 'rms', so update man page
  references.


# gmodels Version 2.15.1 - 2011-01-16

Bug fixes:

- Fix warnings reported by new versions of R CMD check.


# gmodels Version 2.15.0

New features:

- Add support for 'mer' model from lme4.

Bug fixes:

- Correct several minor .Rd syntax errors

- Move extra copyright text to Author field instead of License field.


# gmodels Version 2.14.1

New features:

- Add support for 'lme' objects to estimable().

Other:

- Fix minor typos in manual page for estimable().

# gmodels Version 2.14.0

New Features:

- Add support for 'mlm' objects to estimable

# gmodels Version 2.13.2

Bug Fixes:

- Lower and upper end of confidence interval for lmer objects were
  reversed in est.lmer().

- Correct Greg's email address in two help files.


# gmodels Version 2.13.1

Bug Fixes:

- Problem: R CMD check errors under development version of R 2.5.0
  Solution:
	- Add additional packages to 'Suggests' list in DESCRIPTION
	- Remove extra trailing comma in function calls
	- fix various code/doc inconsistencies

- Problem: estimable() was failing for lmer objects.
  Solution:
	- Create a generic estimable()
	- Move old function to estimable.default()
	- Add  estimable.lmer() to the exported methods list in NAMESPACE

# gmodels Version 2.12.0

- Updated Greg's email address.

- Add support for lmer (lme version 4) objects to ci(), estimable(),
  and fit.contrast() via code contributed by Randall C Johnson.

- Add simplfied coefficient specification to estimable() based on a
  function provided by Randall C Johnson.  It is now possible to do
  things like:
	estimable(reg, c("xB"=1,"xD"=-1))
  instead of:
        estimable(reg, c(    0,   1,	 0,   -1))
  which should make estimable() much easier to use for large models.

# gmodels Version 2.0.8

 - Added DESCRIPTION and removed DESCRIPTION.in

 - Updated CrossTable.R

 - Updated NAMESPACE file

