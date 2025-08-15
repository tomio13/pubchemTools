#/usr/bin/env python
""" Query the pubchem database, and create a record compatible to FURTHRmind
    There are two ways to go: either use the pubchempy module, or go drectly
    pulling the record from the website.
    The module has a nice API, but it is less informative when it comes to
    not so trivial searches.
    Thus, this module adds a bit more to the pubchempy.

    Authors: Tomio
    Date:    2022-10-06-
    License: CC-BY(4)
    Warrany: None
"""
# handle pubchem entries the hard way:
import json
import requests

# some tricks to manipulate the dict tree
from dictDigUtils import dict_search_in_key
#from . ghs_code import *

__all__ = ['Pubchem', 'get_value', 'get_cid']


class Pubchem():
    """ an envelop for pubchem searches
        Use it by searching for an object. Its init will use a search,
        by whatever is requested...
    """

    def __init__(self, search_string: str ='', search_type: str ='name'):
        """ parameters are passed to search direclty to search
        """
        self._record_ = {}

        if search_string:
            self._record_= search(search_string, search_type, 'all')
    # end __init__()


    @property
    def to_dict(self)->dict:
        """ convert all class values to a dict structure,
            with names as keys and values as values.
        """
        drop_list = ['CAS', 'to_dict']
        names = [i for i in dir(self) if i[0] != '_' and i not in drop_list]
        # here we may also employ
        # names.sort()
        # but it is not that important...
        # now, make a dict:
        return {i: getattr(self, i) for i in names}


    @property
    def molecular_weight(self)->float:
        """ Molecular Weight
        """
        res = get_value('Molecular Weight', self._record_)
        if res:
            return float(res[0])

        return -1.0


    @property
    def ghs(self)->list:
        """ all safety codes (H and F codes)
            accessed by searching for the actua H..., P...
            strings in the values of the data tree.
        """
        return get_ghs(self._record_)

    @property
    def translage_ghs(self)->list:
        """ translate the codes to meaningful text
            based on information downloaded from GHS at
            https://pubchem.ncbi.nlm.nih.gov/ghs/
        """
        indx = self.ghs
        res= dict()
        keys=['H-codes', 'P-codes']
        for i in keys:
            if i in indx and indx[i]:
                res[i] = [ghs[k] for k in indx[i]]
        return res

    @property
    def cas(self)->list:
        """ chemical abstracts service codes as list
            Searching for CAS, Related CAS and Deprecated CAS
            in the whole record.
        """
        return get_value_filtered('CAS', self._record_, ['Related CAS', 'Deprecated CAS'])


    @property
    def density(self)->list:
        """ Get the density out of the data.
            Unfortunately, there are many variant of this,
            it may be provided as:
            - a number
            - a number and unit
            - < a value
            - reference to a table image
            - a text describing value ranges and conditions

            For first, we just return the list we find,
            and filter for vapor density
        """
        return get_value_filtered('Density', self._record_, ['Vapor'])


    @property
    def name(self)->str:
        """ RecordTitle from the pubchem record
        """
        if 'RecordTitle' in self._record_:
            return self._record_['RecordTitle']

        return ''


    @property
    def cid(self)->int:
        """ pubchem record ID, cid
            Actually the RecordNumber.
        """
        if 'RecordNumber' in self._record_:
            return self._record_['RecordNumber']

        return -1


    @property
    def iupac_name(self)->str:
        """ IUPAC name
        """
        res= get_value('IUPAC', self._record_)
        if res:
            return res[0]

        return ''


    @property
    def inchi(self)->list:
        """ The InChI of the chemical
            A list, where one field is the InChiKey, the other is
            the actua InChI
        """
        res = get_value('InChI', self._record_)
        if res:
            return res

        return ''


    @property
    def molecular_formula(self)->list:
        """ Molecular formula, as list, because there may be variants
        """
        return get_value('Molecular Formula', self._record_)


    @property
    def smiles(self)->str:
        """ The SMILES descriptor of the molecule
        """
        res= get_value('SMILES', self._record_)

        if res:
            return res[-1]

        return ''


    @property
    def synonyms(self)->list:
        """ the list of synonyms from the record
        """
        return get_value('Synonym', self._record_)
# end class Pubchem


def search(search_string: str ='', search_type: str ='name', what: str ='cids')->dict:
    """ Run a search in the pubchem database based on direct web call.
        Craft an URL, and get the results.
        Return information based on what. If 'cids', then a list of CIDs,
        if 'all', then a full record.
        The full record is cleaned up a bit for easier access and search.
        Use the helpers to find out actual values...

        About synonyms: all are found as name, but some are tricky.
        For example Salzsäure is written as Salzsaeure...

        Parameters:
        search_string:  information to search for, e.g. chemical name
        search_type:    type of information, e.g. name, CAS, formula, cid
        what:           what to return? E.g. cids or record or all (all details)

        Return:
        list of found Pumbed IDs = CID values
    """
    if search_string == '':
        print('Empty query')
        return []

    # we get the full record using the pug-rest API
    domain = 'compound'

    mainlink = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/{domain}/{search_type}'
    # print('link:', mainlink)

    result = {}
    # we want to treat different if what is 'all'
    # so we use an indicator, and we do a search for 'cid's
    # we can use to retrieve the data or send to FURTHRmind
    all_what = False

    if what == 'all':
        all_what = True
        what = 'cids'

    with requests.request(method='GET',
                        url=f"{mainlink}/{search_string}/{what}/JSON",
                          timeout= 30) as a:
        # We have sent the request, what did the server reply?
        # status_code == 200 --> all fine, we got a meaningful content back
        # all others are some kind of reasons why we did not ...
        if a.status_code == 200:
            result= json.loads(a.text)
        else:
            print('call returned:', a.status_code, a.text)
    # end calling the API for CIDs

    if (what == 'cids'
        and not all_what
        and 'IdentifierList' in result
        and 'CID' in result['IdentifierList']):

        return result['IdentifierList']['CID']

    if all_what:
        # refine the search to a more complete one
        # but now we can ride another API for the full record
        # print(result.keys())

        if result == {}:
            print("Nothing found")
            return result

        res = result['IdentifierList']['CID']

        if len(res) > 1:
            print('we restrict to the first full record')

        url = "https://pubchem.ncbi.nlm.nih.gov/rest/"\
              f"pug_view/data/compound/{res[0]}/JSON/?"\
              "response_type=display"

        with requests.request(method= 'GET', url= url, timeout= 30) as a:
            result= json.loads(a.text)

        if a.status_code == 200 and 'Record' in result:
            # result = dict_flatten(clean_section(result['Record']))
            result = clean_section(result['Record'])
        else:
            print('For the full request, server responded:', a.status_code)
            result = {}

    return result
# end search


def get_cid(search_string):
    """ a shortcut to search, returning the CID of a compound based on
        name or CAS number.

        Paramter:
        search_string: text to be searched for, name or CAS number typically

        Return:
        a list of found Pubchem CIDs
    """
    return search(search_string, 'name', 'cids')
# end get_cid


def dig_value(info)->dict:
    """ pubchem records tend to have a 'Value' key,
        under which a list of dicts list up values
        with details for references, and value types
        like String, Numeric, Boolean...
        Dig into these to get the actual values out, because
        python allow to have any of these without extra specifications.
        This makes the structure of the record simpler.

        parameters:
        info:       key within a dict of dicts

        return:
        a dict {'value': value}
    """
    res = {}


    if isinstance(info, str):
        res['value'] =  info

    if isinstance(info, dict):
        if ('StringWithMarkup' in info
            and 'String' in info['StringWithMarkup'][0]):

            res['value']= [i['String'] for i in info['StringWithMarkup']]

        elif 'Number' in info:
            if (isinstance(info['Number'], list)
                and len(info['Number']) == 1):

                res['value']= info['Number'][0]
            else:
                res['value'] = info['Number']

        elif 'Boolean' in info:
            if (isinstance(info['Boolean'], list)
                and len(info['Boolean']) == 1):
                res['value'] = info['Boolean'][0]
            else:
                res['value']=info['Boolean']
        else:
            # print('Unknown value structure!')
            return info

    if 'value' in res:
        return res
    # else
    return None
# end dig_value


def pop_dict_key(info: dict)->str:
    """ Take a dict and look for a specific key in it,
        pop and return this key.
        Potential keys:
        TOCHeading, Name, ReferenceNumber

        This helper function is used to get a key for
        list of dicts, like a section, to turn the
        list to a dict with the key we provide here

        @parameter info:    the element in a section

        @return: the key found
    """
    pop_list = ['TOCHeading', 'Name', 'ReferenceNumber']
    for i in pop_list:
        if i in info:
            return i

    return ''
# end pop_dict_key


def clean_section(info)->dict:
    """ Take a dict returned by Pubchem (within the Record field) and scan it,
        extract the Section lists and put them into the original dict with keys
        obtained from the TOCHeading fields.
    """
    update_list = ['Section', 'Information', 'Value']

    res = {}
    if isinstance(info, dict):
        for k, v in info.items():
            if k.lower() == 'section':
                res.update(clean_section(v))
            else:
                res[k] = v

    elif isinstance(info, list):
        # a list should be a list of dicts

        for j,i in enumerate(info):
            # what shall be a key?
            # to add the list elements to the root dict,
            # we need a key... Candidates are in the pop_list
            # Here we cannot deal with list elements that are
            # not dicts!
            if not isinstance(i, dict):
                continue

            # we hunt for a specific key, its value
            # is used as key for the whole element, and be dropped
            k = pop_dict_key(i)
            if k:
                k = i.pop(k)
            else:
                k = str(j)

            # update what to be updated
            for update_i in update_list:
                if update_i in i:
                    i_subdict = i.pop(update_i)
                    new_value = clean_section(i_subdict) if update_i!='Value'\
                            else dig_value(i_subdict)

                    i.update(new_value)
                    break

            # now, store the result:
            res[k] = i
    else:
        print("Unknown data", info)
    return res
# end of clean_section


def get_value(search_text: str, data: dict)->list:
    """ search for search_text in data, using dict_search_key,
        and then dig for value fields within and return their
        content.

        @param search_text: the key to look for
        @param  data:       the pubchem_search data to dig into

        @return: a list of hits
    """

    res_list = dict_search_in_key(search_text, data)

    res = []
    for v in res_list:
        if isinstance(v, dict):
            res_line = dict_search_in_key('value', v)

            if res_line:
                if isinstance(res_line[0], list):
                    res += [i for j in res_line for i in j]
                else:
                    res += res_line

    # make the results unique
    if res:
        return list(set(res))

    # else
    return []
# end get_value


def get_value_filtered(search_text: str,
                       data: dict,
                       kill_list: list) -> list:
    """ Use the get_value above but filter the resulted
        keys for ones in kill_list, and drop those listed
        there.


        @param search_text: the key to look for
        @param  data:       the pubchem_search data to dig into

        @return: a list of hits
    """
    if not kill_list:
        return get_value(search_text, data)

    kill = []
    for i in kill_list:
        kill += get_value(i, data)

    res = get_value(search_text, data)

    if kill:
        for i in kill:
            if i in res:
                res.remove(i)
    return res
# end get_value_filtered


def get_ghs(info):
    """ Dig out the GHS codes, that is the 'H' and 'P' codes.
        This is a bit brute force analyzing the values, that is
        the descriptions themselves, not the actual dict structure.

        @parameter info the pubchem record
        @return: dict with 'H-values' and 'P-values'
    """

    if not isinstance(info, dict):
        raise ValueError('A pubchem record dict was expected')

    codes = get_value('GHS', info)
    # this list should have two types of values:
    # Hxxx: explanation
    # Pxxx, Pxxx, Pxxx.... list
    # some numbers occasionally, like record reference

    h_codes = []
    p_codes = []
    for i in codes:
        if isinstance(i, str):
            if i.startswith('H') and ':' in i:
                h_codes.append(i.split(':')[0])
            elif i.startswith('P') and ',' in i:
                p_codes += [j.strip() for j in i.split(',')]


    # last P value starts with 'and ...'
    if p_codes and  p_codes[-1].startswith('and'):
        p_codes[-1] = p_codes[-1].split('and ')[-1]

    return {'H-codes': h_codes, 'P-codes': p_codes}
# end of get_ghs


# you can test running:
# c = pubchem_search('sulphuric acid', 'name', 'all')
# cas = get_cas(c)
# id = c['RecordNumber']
# m = get_value('Molecular Weight', c)[0]
# form = get_value('Molecular Formula', c)[0]
# Hazard codes:
# codes = get_ghs(c)
