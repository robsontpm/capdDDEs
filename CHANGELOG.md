
Change Log
==========

All notable changes to this project will be documented in this file,
especially the ones pertaining to the ***rigorous*** nature of
the computations. If any serious bug is found affecting that, it
will be reported here.  
 
The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project adheres to [Semantic Versioning](http://semver.org/).

 
## [Unreleased] - 2024-12-13
 
### Fixed

 - (serious) a bug that was caused in CAPD with the AutoDiff 
   fadbad library and capd::vectalg::Vector the bug was introduced
   after using the submodule CAPD instead of build-in version. 
   Only codes using version from 2017-11-15 were affected. 
 - a bug introduced during the process of debugging the aforementioned
   autodiff bug. 
 - (serious) a bug in solving non-autonomous systems. If you have ever
   used this library to solve non-autonomous DDEs, please 
   check your codes. 

 
## [Unreleased] - 2017-11-15
 
### Changed

Moved from internal capy of CAPD to a submodule of newest 
CAPD from the official github repository. Please follow
the steps described in READM.md to update the submodule. 
The project is now public and everybody should use this 
version instead of codes taken from my webpage. 

 
## [Unreleased] - 2017-02-19
 
### Added

Initial commit, copy of the old codes available at my 
webpage [scirsc.org](http://scirsc.org). repo is private for now. 