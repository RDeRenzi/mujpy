.. _philostory:

Philosophy & History of MuSR code
==========================================
Just a bit. This is a place for some speculations. 

The killer of data analysis software is sustainability (i.e. who maintains it)
In MuSR  very good software packages, like <wimda>_, are mostly a one man effort.  <musrfit>_ may be an exception, but it is certainly not an open sofwtare effort. Both have the great advantage of a very professional user as the maintainer (the guru), therefore what you need is there, although through the personal preference of the guru. 

Many attempts to build a software as a professional collective fight principally with the reduced community of user, not justifying very large efforts, but also with professional programmers lacks of the *feeling* of real needs.

Some of the needs are a bit complex. The user should find it extremely easy to think of the model and to modify it to suit the physics under investigation.

A crucial part of this is global fits. They pose two problems: they are cumbersom to build and adjust; their goodness is not straightforward to assess. <musrfit>_ offers a very good solution, but admittedly it has a *steep learning curve*.

Another issue is ease of installation. Certainly CERN root is the source of minuit, hence a necessary reference for physics data fits. But it is a bit of a *bloatware*.


In this context herw we explore another solution seeking the following advantages

* Intuitive GUI approach, with attention to power-execution, presently from single fits to suites.
* Based on python, a present growing standard.
* Inserted in Jupyterlab, a growing standard for data analysis, and offering with that a Browser interface
* Easy to install, also on macOS and Windows.
* Easily documented by ReadTheDocs
* Implementing continuous integration [not yet ;)] 
  
