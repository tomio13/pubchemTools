# pubchem tools
is a small packet to query the pubchem APIs. For simple and easy access,
one can use the pubchempy library.
It is not as nice as pubchempy neither so complete, but works well for
things like finding CAS numbers, safety codes, and anything one knows in
the record.
Key point: it works with pulling the record, not querying it every time
one asks a question.

## pubchem has two interfaces
Actually the PUG-RES API is available two ways:
 * [to query a compound by name](https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/)
 * [or to get a full record in JSON format](https://pubchem.ncbi.nlm.nih.gov/rest/pug\_view/data/compound/{res[0]}/JSON/?response\_type=display)

and they provide slightly different facilities.
Using the former, one can get simple answer to things like:
what is the CID (or pubchem ID) of compound 'a'?
The latter fetches the full record in a hierarchic JSON object.

## using only pug
Trying the various output methods did not provide a response as rich as
the pug\_view did. Thus, the double call is important at the moment.

# Handling the JSON hierarchy
One thing to do is to pick out the internal lists of dicts to
some objects, e.g. turning them to dicts. Then a crawler
across the dict structure can discover information we need.
We use the dictDigUtils, a c.a. 5 kB set of functions help digging the
hierarchy returned by the json python library.

Lists are:
* Section
* Value
* StringWithMarkup
* Number
* Boolean

Funcitons like clean\_section help managing this.

A Pubchem class is formed which has simple get\_value calls to provide
parameters, such as molecular\_weight, molecular\_formula, etc.

# full record
is available as the \_record\_ private variable. Some parts are simplified
during the processing to remove formating instructions, not needed
here.
You can always use:

```python
a = Pubchem('your molecule')
get_value('Molecular Weight', a._record_)
```

to query into the record. dictDigUtils can provide you dict\_list\_keys
to discover possible keys within.

# names and synonyms
The PUG-RES API will recognize any names what is listed in the synonym list,
thus one needs no other trics to find the material. However, foreign names
may be encoded e.g. the German ä as ae.
Now the list of synonyms is available as the synonyms element of the class.

# InChi
There are two fields coming from Pubchem:
 * InChI key
 * InChI, which is indicated with InChi= in the text

# density
The data set may contain several versions of this, and unfortunately
it is not really uniform across the chemicals.
There are:
 - a number
 - a number and unit
 - < a value
 - reference to a table image
 - a text describing value ranges and conditions
So, for now, the function just returns the list it finds.

A bit more generic solution is to allow a filter in the general searc.
This is now the get_value_filtered() function, which handles the density
and the CAS numbers. For density, use a cut for vapor values. Others did
not work, because they were listed under the same key.

## 1.0.9
Fix a typo in pubchem CAS field, coming from the update for density.

### add list of H and P code meaning
downloaded from https://pubchem.ncbi.nlm.nih.gov/ghs/, and extracted only
the textual meaning part for better human readability.
And changed the license to CC-BY(4)
Update the dependency list in setup.cfg
Update the docstrings to be a bit more informative.

