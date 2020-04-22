# Pretty Good Software Engineer

![Little red dougie hood](https://i.imgur.com/GXTiRn4.jpg)

## Source code for for my personal website and blog

This website is hosted using [Github Pages](https://pages.github.com/) at https://graves.github.io. The blog is built using the wonderful [mdBook](https://github.com/rust-lang/mdBook).

> mdBook is a utility to create modern online books from Markdown files.

The index page and individual posts are written in [Markdown](https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet), compiled using mdBook, and it's output pushed to the [graves.github.io repository](https://github.com/graves/graves.github.io).

The result is a beautiful blog in the style of Rust's documentation.

![Screenshot of blog](https://i.imgur.com/QiMtzdC.png)

Additional language highlighting for codeblocks (Smalltalk in particular) was added by compiling my own [highlight.js](https://highlightjs.org/) file and replacing the existing file in the theme directory.

## Hacking / Contributing

To begin hacking on the blog (to perhaps use as a starting point for your own) you must first [install Rust](https://www.rust-lang.org/tools/install) and then mdBook.

`cargo install mdbook`

Once the prerequisitees are out of the way hacking is as simple as cloning the repository

`git clone https://github.com/graves/pretty_good_software_engineer`

 and adding your own Markdown files to the [src](src) directory.

`head -n 7 src/covid-19-SIR-model.md `
```markdown
# Modeling the Transmission Rate of COVID-19 with SIR in Smalltalk

## Using Glamorous Toolkit, PolyMath, Roassal2 and Ordinary Differential Equations

From [Wikipedia](https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology#The_SIR_model):
> The SIR model is one of the simplest compartmental models, and many models are derivatives of this basic form. The model consists of three compartments: S for the number of susceptible, I for the number of infectious, and R for the number of recovered or deceased (or immune) individuals. This model is reasonably predictive[citation needed] for infectious diseases that are transmitted from human to human, and where recovery confers lasting resistance, such as measles, mumps and rubella.
```
 
 The posts can then be referenced in [src/SUMMARY.md](src/SUMMARY.md)

`cat src/SUMMARY.md`
 ```markdown
# Summary

[Home](./home.md)
- [Modeling the Transmission of Covid-19 with Smalltalk](./covid-19-SIR-model.md)
```

Now you can build the website

`mdbook build`

The result with be output into a directory named `book`. The contents of which can then be uploaded to and served by your favorite HTTP server or served locally

`mdbook serve`

You can then navigate to your site at http://localhost:3000
