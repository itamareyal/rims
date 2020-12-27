'''
ion_parser.py
'''


import periodictable

def ion_inputs_to_attributes( str_ion):
        ion_lst = []
        ion_dict = []

        if ',' in str_ion:
            splitted = str_ion.split(',')
            for ion in splitted:
                ion_lst.append(ion)
        else:
            ion_lst.append(str_ion)

        for ion in ion_lst:
            ion = ion.strip(' ')
            splitted_ion = ion.split(' ')
            symbol = splitted_ion[0]
            # check that the symbol is legit
            if splitted_ion[1][0] in ['+','-']:
                charge = int(splitted_ion[1][1:])
            else:
                charge = int(splitted_ion[1])
            mols = float(splitted_ion[2])

            if splitted_ion[1][0] == '-':
                charge *= -1

            ion_dict.append({
                "symbol": symbol,
                "charge": charge,
                "mols": mols
                })
        return ion_dict


def verify_ions(ions_dict_lst):
    for el in ions_dict_lst:
        symbol = el["symbol"]
        try:
            element = periodictable.elements.symbol(symbol)
        except:
            print(symbol +" wasn't found in periodictable lib.")
            return False
        el["element"] = element
    return True


def print_ions(ions_dict_lst):
    print("\nThe following ions were added to the system:")
    for el in ions_dict_lst:

        if el["charge"] > 0:
            sign = '+'
        else:
            sign = ''

        print(el["element"].symbol +" " + el["element"].name+"\n\tcharge: "+sign+str(el["charge"])+"[Columb]\n\tquantity: "+str(el["mols"])+" [mol]\n")

    input("\npress ENTER to confirm...\n")



