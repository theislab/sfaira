---
name: Data loader request
about: Suggest a new dataset to add to sfaira.
title: 'Dataset summary'
labels: 'dataset'
assignees: ''

---

**Is your dataset not yet in sfaira?**
Please search the sfaira [code](https://github.com/theislab/sfaira/tree/dev), and [issues and PRs](https://github.com/theislab/sfaira/issues?q=) for any DOI associated with this dataset, such as preprints or journal publications. To be safe, you can also search for the first author.

**Describe the data set**
Fill in below:

- Link to publication:
- Title of publication:
- First author of publication:
- DOI preprint:
- DOI journal:
- Download link to any data objects required:

**Describe meta data**
Optionally, you can already collect meta data information here before starting to write the data loader, or to help others write one.
You can extend these lists if you find more meta data that you want to record before writing a data loader.
Note that you can also directly put this information into a draft data loader and start a pull request instead of first writing an issue.
If you know this dataset well but cannot write a loader right now, this will help other people decide if they want to write a loader for this dataset and will speed up their process.

- Is this primary data (not a meta study):
- Is most raw gene expression matrix normalized (if yes how)?:
- Single-cell assay used:
- Disease(s) of sampled individuals (ideally MONDO term):
- Organ(s) sampled (ideally UBERON term):
- Organism(s) sampled (ideally NCBItaxon term):

Any relevant cell-wise annotation fields that are column names of a table or column names in `.obs` of an `h5ad` for example:

- Cell type annotation:

**Additional context**
Add any other context of the dataset or anticipated issues with the dataloader.
