#!/usr/bin/env python
#import configs

class configs():
    '''Read configuraitons and do stuff.

        Input::train_x = file contains configuraitons in format.
                1. Sample:

                    6
                       -0.00476416    -0.90668842     0.00000000     0.00000000
                    O   1.39511118877395       -1.15984255934314       4.127641665696761E-015
                    O   1.39511118877394        1.15984255934314      -1.010148437748930E-014
                    O  -1.37600000000000      -4.445978668838478E-015  2.222989334419239E-015
                    H  -1.99243331768939      -0.734636620951108       6.321962738314389E-015
                    H  -1.99243331768940       0.734636620951091       1.574609909686252E-015
                    C   1.37600000000000       0.000000000000000E+000  2.222989334419239E-015

                2. Mapping to parameters:




                    Breakdown:
                        configs =

                            [
                                [config],    #configs[0], first configuraiton
                                [config],    #configs[1], second configuraiton
                                ...
                                [config]     #configs[-1], last configuraiton
                            ]


                        config = configs[0] =    #Take first configuraiton for example
                            [
                                [molecule_count_total],    #number of atoms in a molecule
                                [[energy], [dipole]]
                                [molecule]
                            ]
                                    #molecule_count_total (int):    configs[0][0] = 6
                                    #energy (float)(Hartree):    configs[0][1][0] = -0.00476416
                                    #dipole (float)(a.u.): configs[0][1][1][0]=-0.90668842


                        molecue = configs[0][2] = config[2] =
                            [
                                [[element],   [x,   y,  z]],   #First atom
                                [[element],   [x,   y,  z]],   #Second atom
                                ...
                                [[element],   [x,   y,  z]]    #Last atom
                            ]
                                    #element  (string) : configs[0][2][0][0]
                                    #x/y/z (float)(Angstrom):   configs[0][2][0][1][0] = '1.39511118877395'

    '''

    def __init__(self,train_x,dip = False,first_n_configs=False):
        '''Read input file into configs with checks

        1. Result:  configs[a][b][c][d][e], type = List
            configs stores the all the input obtained from train_x
                User's guide:
                [a]:  The ([a]+1)_th configuration #Python list count from 0

                    [b]=0: Number of atoms

                    [b]=1:

                        [c]=0: energy
                        [c]=1: dipole

                            [d]=0: dipole_x
                            [d]=1: dipole_y
                            [d]=2: dipole_z

                    [b]=2: molecule

                        [c]: The ([c]+1)_th atom

                            [d]=0: element
                            [d]=1: cartisian coordinate

                                [e]=0: x
                                [e]=1: y
                                [e]=2: z

        2. Check input:
            While reading files, it can check:
            * Skip blank lines automatically.
            * Number of atoms: existence, type and consistency
            * Energy: existence and type
            * Dipole (if dip = True): existence, type and dimension
            * Element: existence, type and consistency
            * Cartisian coordiante: existence, type and dimension
            Notice: consistency error includes the existence error because not existent is also considered as inconsistent.
        '''
        """Constants"""
        import copy
        self.hartree_to_cm = 219474.63


        self.logo()
        print('Reading file...\n')
        f = open(train_x)
        configs_count = 0
        line_count = 0
        element_1 = list()
        configs = list()
        dipole = list()
        blank_line_count = 0
        for line in f:
            line = line.strip() #The delete the newline character
            #if len(line) == 0: #Comment out so that don't can if E is 0
               # blank_line_count += 1
                #continue
            line_count = line_count + 1
            if line_count == 1: # This is number of atoms.  #Check type
                if line.isdigit():
                    molecule_count_total = int(line)
                    #break
                    config_count_totol = molecule_count_total + 2
                else:
                    print("Type error: molecule number in line "+str(line_count)+" is : "+line)#Check type

            config_count = (line_count - 1) % config_count_totol + 1

            if configs_count == 1 and config_count >= 3 and config_count <= molecule_count_total+2:
                element_1.append(line.split()[0]) #Read elements in first config

            if config_count == 1:
                molecule_coord = list() #Renew molecue_coord each configuration
                configs_count = configs_count + 1
                try:
                    if int(line) != molecule_count_total:
                        print('Consistency error: number of atoms in line: ' + str(line_count)+' is : ' + line) #Check consistency
                except:
                    print("Type error: number of atoms in line "+line_count+" is : "+line) #Check input
                    #break

            if config_count == 2:
                try:
                    if len(line) is 0:
                        energy = [0]
                    else:
                        energy = [float(line.split()[0])] #Record the energy.
                except:

                    print("Type error: energy in line "+line_count+" is : "+line) #Check energy type
                    #break


                try:
                    if dip == True and len(dipole) !=3:
                        dipole = [float(dip) for dip in line.split()[1:]] #Record dipole.
                        print("Dimension error: Dipole dimension in line "+str(line_count)+" is : "+str(len(dipole)))
                except:
                    print("Type error: Dipole in line "+str(line_count)+" is: "+str(line)) #Check energy type
                    #break


            if config_count >= 3 and config_count <= molecule_count_total+2:
                molecule_count = config_count - 2
                element = line.split()[0]

                if element != element_1[molecule_count-1] and len(element) != 0: # Element check type and existence
                    print('Consistency error: atom in line: ' + str(line_count)+' is: '+line)

                coordinate = [float(coor) for coor in line.split()[1:]]#Read cooridnates and turns into float
                atom_coord = [element,coordinate]
                molecule_coord.append(atom_coord)

            if config_count == config_count_totol: #Reset some of the values
                configs.append([[molecule_count_total],[energy, dipole],molecule_coord])


                if configs_count == first_n_configs: #Only read first n configs
                    break


        self.blank_line_count = blank_line_count
        self.configs = copy.deepcopy(configs)
        self.dip = dip
        self.configs_count = configs_count
        self.train_x = train_x
        self.molecule_count_total = molecule_count_total
        self.line_count = line_count


        self.configs_sorted = self.sort(configs)
        self.energy_lowest=self.energy_array_sorted[0]
        self.energy_highest=self.energy_array_sorted[-1]
        self.energy_lowest_cm = self.energy_array_sorted_cm[0]
        self.energy_highest_cm = self.energy_array_sorted_cm[-1]
        print('Number of blank lines in file:    {:<2d}'.format(self.blank_line_count))
        print('Configuration check finished!')
        #print('**Status: File reading finshed.\n')
        print('Number of configurations:    {:<6d}'.format(self.configs_count))

        print('Reading finished. You can use self.info() or self.help() to start')
        #self.info()
        #self.help()


        f.close()


    def todo(self):
        print("""
        TODO:

        1. Add feature so that input can be formated string.
        2. Add molden path
                """)

        #    def read(self,config_string):
        #for line in string:

    def logo(self):
        print("""

        ########  ##    ## ########  ########  ######
        ##     ##  ##  ##  ##     ## ##       ##    ##
        ##     ##   ####   ##     ## ##       ##
        ########     ##    ########  ######    ######
        ##           ##    ##        ##             ##
        ##           ##    ##        ##       ##    ##
        ##           ##    ##        ########  ######

                                            Version 0.0.36

                                --A Bowman Group Product
                                    """
                )
        """--Developed by Kee"""

    def info(self,configs=False):


        print('============================')
        print('||          INFO          ||')
        print('============================')

        print('Inputfile:  {}'.format(self.train_x))
        print('Number of atoms:    {:<2d}'.format(self.molecule_count_total))
        self.order()
        print('Number of configurations:    {:<6d}'.format(self.configs_count))

        print('\nLowest energy (Hartree):    {:14.8f}'.format(self.energy_lowest))
        print('Lowest energy (wavenumber):    {:14.8f}'.format(self.energy_lowest_cm))
        self.prt(self.configs_sorted[0])
        print('Highest energy (Hartree):   {:14.8f}'.format(self.energy_highest))
        print('Highest energy (wavenumber):   {:14.8f}'.format(self.energy_highest_cm))
        print('Configuration that contains highest energy: ')
        self.prt(self.configs_sorted[-1])
        print('========END OF INFO=======\n')

    def feature(self):
        print('''
            Those are important features:
                1. Sorted configs according to energy and pretty writing.
                    a=configs('input.xyz')
                    a.write('output.xyz',a.sort())

                2. Sorted configs according to energy and print first 10 on screen.
                    a=configs('input.xyz')
                    b=a.list()[0:10] #Turn all a configs into a list and choose first 10 of it
                    a.prt(b) #Pretty printing

                3. Sort configs and see then using molden.
                    a=configs('input.xyz')
                    a.write('output.xyz',a.sort())

                4. Resize distance between monomers.
                    a=configs('input.xyz')
                    a.resize() #(Then follow the structure)

                *5. Plot energy distribution with binwith (defaut is 50 cm-1).
                    a=configs('input.xyz')
                    a = configs('input.xyz')
                    a.plot(binwith=10)

                *Might not work well. Still testing.

                ''')

    def order(self,monomer = False, configs=False):
        """TODO: Add feature to show the group of monomers"""


        configs=self.configs_check(configs)
        print('Atom numbering:')
        molecule_count = 1
        for atom in configs[0][2]:
            print('{} ({})'.format(atom[0],molecule_count))
            molecule_count = molecule_count + 1

    def help(self):
        syntax="""
        =============================
        ||       SYNTAX HELP      ||
        =============================

        Suppose you have:
            a = configs('input.xyz')

        Now you can:

            a.prt(configs=False)
                *Print the configs to the screen with pretty format

            a.write('output.xyz', configs=False)
                *Write the list of configs into output file with pretty format

            a.plot(self,configs = False,binwidth=False)
                *Plot energy distribution with gien binwith (default is 50 cm-1)

            a.sort(configs=False,reverse=False,key ='energy')
                *Sort the configs according to energy and return the list

            a.switch()
                *Switch the order to atoms of all configs as return as list.

            a.list(first_n_configs=False)
                *Return default configurations as a list

            a.threshold_energy(configs,lower=False,upper=False)
                *Select all configurations that satisfy the energy threshold

            a.distance(config,atom_A,atom_B)
                *With given one configuraiton, calculate the distance of two atoms

            a.order()
                *Show the order of atoms

            a.resize(self, n_configs=False,configs=False, first_n_configs=False, monomer=False, dis_lower=False, dis_uppder=False, dis_new_lower=False, dis_new_upper=False)
                *Select configs randomly in the given region[dis_lower, dis_upper], then separate them along given two atoms to the given new distance [dis_new_lower, dis_new_upper].

            a.dissociation(self,config=False, dis_min=False,dis_max=False,step=False, atomA=False,atomB=False)
                *Dissociate two monomers with given

        *To show info:
            a.info()

        *To show syntax help:
            a.help()

        =======END OF HELP SYNTAX====

        """
        #for line in syntax.split('\n'):
        #            print(line.strip())
        print(syntax)

    def prt(self,configs=False):
        """Print configs into screen.

        """

        #print('*Status: printing...')
        configs_count=0
        configs=self.configs_check(configs)

        for config in configs: #Write config by config
            configs_count = configs_count + 1
            print('    {:<2d}'.format(config[0][0])) #Number of atoms. Align number of atom to the very left
            print(' '),# To align the colums of energy and coordiante
            if self.dip:#Print option for having dipole or not
                print('    {:14.8f}{:14.8f}{:14.8f}{:14.8f}'.format(config[1][0][0],config[1][1][0],config[1][1][1],config[1][1][2]))
            else:#Print also dipole (if input file has it)
                print('    {:14.8f}'.format(config[1][0][0]))#Print energy only
            for molecule in config[2]:#The atom part: coordiante of atoms
                print('    {} {:14.8f}{:14.8f}{:14.8f}'.format(molecule[0],molecule[1][0],molecule[1][1],molecule[1][2]))

        print('----Printed {} configuration(s).\n'.format(configs_count))
        #print('*Status: printing finished.\n')

    def write(self,filename, configs=False): #Default value has to be after non-default value.
        """Write formated configurations into filename.

        """

        print('*Status: writting file...\n')
        f=open(filename,'w')
        configs_count=0

        configs=self.configs_check(configs)
        for config in configs: #Write config by config
            configs_count = configs_count + 1

            f.write('{:<2d}'.format(config[0][0])+'\n') #Number of atoms. Align number of atom to the very left
            f.write('  '), #To align the colums of energy and coordiante. Comma to continue on the same line
            if self.dip:#Print energy and dipole
                f.write('{:14.8f}{:14.8f}{:14.8f}{:14.8f}'.format(config[1][0][0],config[1][1][0],config[1][1][1],config[1][1][2])+'\n')
            else:#Print energy only
                f.write('{:14.8f}'.format(config[1][0][0])+'\n')#Print energy only
            for molecule in config[2]:#The atom part: coordiante of atoms
                f.write('{} {:14.8f}{:14.8f}{:14.8f}'.format(molecule[0],molecule[1][0],molecule[1][1],molecule[1][2])+'\n')
        print('----{} configs are written to the file: {}\n'.format(configs_count,filename))
        print('*Status: writting finished.\n')
        f.close()

    def sort(self,configs=False,reverse=False,key = 'energy',subkey1=None,subkey2=None): #sort based on energy
        """Sort given configurations based on the given key and return the list.

        TODO: sort key = distance
        """
        import numpy as np
        import copy
        print('*Status: sorting...\n')
        configs=self.configs_check(configs)

        if key == 'energy':
            print('\n----Sorting by energy...\n')
            if reverse: #From high to low
                configs.sort(key= lambda item:item[1][0],reverse = True) #In key, it helps iterate though list.
                #print('\nReverse sort finished.\n')
                #return configs
            else: #Default, from low to high.
                configs.sort(key= lambda item:item[1][0])
                energy_array_sorted = list()
                energy_array = list()

                for config in configs:
                    #print(config[1][0][0])
                    energy_array.append(config[1][0][0])
                self.energy_array_sorted = energy_array
                energy_array_sorted_cm = copy.deepcopy(self.energy_array_sorted)
                self.energy_array_sorted_cm = np.array(energy_array_sorted_cm)* self.hartree_to_cm

        elif key == 'distance': #Reserved for other key expansion in the future.
            print('\n----Sorting by distance...\n')
            if reverse: #From high to low
                configs.sort(key= lambda item:item[1][0],reverse = True) #In key, it helps iterate though list.
                #print('\nReverse sort finished.\n')
                #return configs
            else: #Default, from low to high.
                configs.sort(key= lambda item:self.distance(item,subkey1,subkey2))
                energy_array_sorted = list()
                energy_array = list()

                for config in configs:
                    #print(config[1][0][0])
                    energy_array.append(config[1][0][0])
                self.energy_array_sorted = energy_array
                energy_array_sorted_cm = copy.deepcopy(self.energy_array_sorted)
                self.energy_array_sorted_cm = np.array(energy_array_sorted_cm)* self.hartree_to_cm



        print('\n*Status: sorting finished.\n')
        return configs

    def sort_distance(self,configs=False,reverse=False):
        return

    def list(self,first_n_configs=False):
        """Return a list of reading file.

        """


        if first_n_configs == False:
            return self.configs
        else:
            return self.configs[0:int(first_n_configs)]

    def threshold_energy(self,configs,lower=False,upper=False):
        """First sort the list then select the configs based on energy threshold

        """

        print('Excuting energy_threshold:----------------')
        configs = self.sort(configs)
        configs_new_count = 0
        configs_new = list()
        energy_lowest=configs[0][1][0][0]
        energy_highest=configs[-1][1][0][0]


        if lower is not False:
            energy_lowest = float(lower)
        if upper is not False:
            energy_highest = float(upper)
        if energy_lowest > energy_highest:
            print('Boundary error: lower energy bound is greater then upper energy bound. Returning origianl configurations.')
            self.error_message()
            return configs
        else:
            for config in configs:
                if config[1][0][0] >= energy_lowest and config[1][0][0] <= energy_highest:
                    configs_new.append(config)
                    configs_new_count = configs_new_count + 1
            percentage = round(float(configs_new_count)/float(self.configs_count)*100,2)
            print('There are '+str(configs_new_count)+'/'+str(self.configs_count)+'('+str(percentage)+'%) configs between '+str(energy_lowest)+' and '+ str(energy_highest)+' Hartree')
            print('End of energy_threshold-------------------')
            return configs_new

    def threshold_distance(self,configs,lower=False,upper=False):
        return

    def error_message(self):
        """To show error message and waiting for decision
        Internal
        """
        decision = input('Do you want to continue? (y/n) \n')
        if decision is not 'y':
            print('Exiting program.')
            exit()

    def distance(self,config,atom_A,atom_B):
        """To calculate distance with given two atoms.
        Internal
        atom = molecue[i] = [[element],[x, y, z]]
        """

        import numpy as np
        import math
        #print('This is atom_B:',atom_B)
        atom_A = int(atom_A)-1
        atom_B = int(atom_B)-1
        atom_A = np.array(config[2][atom_A][1])
        atom_B = np.array(config[2][atom_B][1])
        dis = math.sqrt(np.sum(np.square(atom_A-atom_B)))
        return dis

    def switch(self,configs=False):
        import copy
        print('Switching order:\nThe origianl numbering of given configs is:')
        self.order(configs)
        molecule_count = 1
        molecule_new_count = 1


        if configs is False:#Default input check module
            configs = self.configs
            print('Switching using orignial configuraitons')

        print('\nWhat new order do you what?\n')

        try:
            temp = configs[0][0][0]
        except:#If has to configuration, this makes sure it can be safely iterated in the next statemnt.
            configs = [configs]

        configs_new = copy.deepcopy(configs)#This is the correct way to create a different list with same value

        while molecule_new_count <= self.molecule_count_total:#Will repeat molecule_count_total times
            configs_new_count = 0

            try:
                molecule_count=int(input('I want new atom ({:d}) to be the old atom number: '.format(molecule_new_count)))
                for config_new in configs_new:
                    configs_new[configs_new_count][2][molecule_new_count - 1] = configs[configs_new_count][2][
                        molecule_count - 1]
                    configs_new_count = configs_new_count + 1
            except:
                pass
            #molecule_count = self.molecule_count_total - molecule_new_count + 1 #Test arguement(reverse order)



            molecule_new_count = molecule_new_count + 1
            #        self.prt(configs_new)
        print('New numbering:')
        self.order(configs_new)
        print('''New configs are returned as list, please use a.write(configs_new,'ouput') to save configs.''')
        return configs_new

    def configs_check(self,configs,silence=True):
        """To check if give the configs arugment.

            Usage: configs=configs_check(configs)
        """

        if configs is False:
            if silence is not True:
                print('Using original list.')
            return self.configs
        else:
            try:
                check = configs[0][0][0]#Check if the list has only one configuration but one layer smaller
                return configs
            except:#If has to configuration, this makes sure it can be safely iterated in the next statemnt.
                return [configs]

    def translate(self,atom_A=False, atom_B=False, dis_new = False, config = False):
        """Give a dimer and designated atom groups, translate in vector AB from previous distance to new distance

            Input:
                dis_new :: the distance after translation, essiential
                monomer :: contains the info about monomer groups and referenced atoms, omit this if use default.
                config :: take one configuraiton so that it can do the translation. If given multiple configuraitons, return only first one.

            Output:
                configs_new :: none-nested configuraiton, subset of configs.

            Theory:
                For point A and B. Now translate B along AB direction to distance dis_new = |AC| and become C. Parameter function for line AB:
                    xC = xB + m*t
                    yC = yB + n*t
                    zC = zB + p*t
                where:
                    t = (dis_new-dis)/sqrt((m**2+n**2+p**2))
                    (m,n,p) = vector(AB)

                Suppose D is also the same monomer as B, and in ordor to translate D:
                    xD_new = xD + m*t
                    yD_new = yD + n*t
                    zD_new = zD + p*t
                where t and m,n,p take the same value as previous one.


        """
        import numpy as np
        import copy
        #print("-------------Atom_A",atom_A)
        configs = self.configs_check(config)
        #print(type(self.monomers()))
        monomers = self.monomers()
        #print('test')
        #print(monomers)

        monomer_A = monomers[0]
        monomer_B = monomers[1]
        #print ('This is monomer_B:',monomer_B)
        if atom_A is False: #Assign atom a and atom b
            try:
                trail = input()
                self.order()
                atom_A = int(input('Please assign the number of first atom: '))
                atom_B = int(input('Please assign the number of second atom: '))
            except:
                print('Translate using default Atom (first atom in monomer_A and first atom in monomer_B)')
                #print(monomer_A)
                atom_A = monomer_A[0]
                #print ('This is atom_A:',atom_A)
                #atom_B = monomer_B[0]
                print ('This is atom_B:',atom_B)
                dis_new = 10



        config_new = copy.deepcopy(configs[0]) #Everytime need a new list with same vaule need deepcopy.
        #print(monomer_A)
        #print('New atom_A: ',atom_A)
        #print('New atom_B: ',atom_B)

        dis = self.distance(config_new,atom_A,atom_B)
        v_AB = np.array(self.vector(config_new, atom_A, atom_B))
        v_A = np.array(self.vector(config_new, atom_A))
        v_B = np.array(self.vector(config_new, atom_B))
        t = (dis_new-dis)/np.sqrt((np.sum(np.square(v_AB))))

        monomer_B_new = list()
        for atom in monomer_B:
            v_atom = np.array(self.vector(config_new, atom))
            v_atom_new = v_atom + v_AB * t
            config_new[2][atom-1][1] = v_atom_new
        return config_new # One none-nested configuraiton list

    def monomers(self, monomer_A = False, monomer_B = False):
        # type: (object, object) -> object
        """To assign monomer groups and return as two lists of number. Use default groups if no assign.

            Example: monomer = [[[3, 4, 5], [6, 1, 2]]]
        """

        monomers=list()

        try:
            if len(self.monomer_AB) is not 0:
                monomers = self.monomer_AB
                 #Meaning monomer has already been assigned. Once assigned forever assigned.
        except:
            print('Monomer is not determined. Please specify monomers')
            #if monomer_A is False: #For the convinience of other function's usage
            try:
                #a = input()#To try out if can catch input
                self.order()
                print('(Enter integers and separate them by whitespace)')
                str1 = input('What atoms are in first monomer: ')
                str2 = input('Waht atoms are in second monomer: ')
            except:
                print('Warning: Using test arguments. Please use terminal to catch input.')
                atom_A= 3# Arguemnt for test
                atom_B= 6# Argument for test
                str1 = '3 4 5'# Argument for test
                str2 = '6 2 1'# Argument for test
            monomer_A = list()
            monomer_B = list()
            for line in str1.split():
                monomer_A.append(int(line.strip()))
            for line in str2.split():
                monomer_B.append(int(line.strip()))
            monomers=[monomer_A, monomer_B]
            self.monomer_AB = monomers
        #print monomers

        #monomer_A = [atom_A]
        #monomer_B = [atom_B]

        return monomers
        #Example: monomers = [[[3, 4, 5], [6, 1, 2]]]

    def dissociation(self,config=False, dis_min=False,dis_max=False,step=False, atomA=False,atomB=False):
        """This method is used to make rigid dissociation along designated atoms.

        """


        if config is False: #Meaning use the default one
            print('Using default global minimum configuration: ')

            config = self.configs_sorted[0] #Default is global minimum
            self.prt(config)

        if step is False:
            try:
                step = float(input('Please specify the step: '))
            except:  #False-save
                print('Using default step: 0.05')
                step = 0.05
        else:
            print('Using default step: 0.05')
            step = 0.05

        if dis_min is False and dis_max is False:
            try:
                dis_min = float(input('Please specify the dis_min: '))
                dis_max = float(input('Please specify the dis_max: '))
            except:
                dis_min = 2
                dis_max = 20

        else:
            print('Using default dis_min: 2')
            print('Using default dis_max: 20')
            dis_min = 2
            dis_max = 20

        configs_count = 0
        monomers = self.monomers(configs)
        (atom_A,atom_B) = self.chosen_atom()
        configs_new = list()
        while dis_min <= dis_max:
            configs_new.append(self.translate(atom_A,atom_B,dis_min,config))
            configs_count = configs_count + 1
            dis_min = dis_min + step
            print('New distance: '+str(dis_min))
        #print(configs_new)
        try:
            decision = input('Do you want to see the configs in Molden? (y/n):')
            if decision is 'y':
                self.molden(configs_new)
        except:
            self.molden(configs_new)


        print('{:d} configs are return as list.'.format(configs_count))
        return configs_new

    def chosen_atom(self,atom_A=False,atom_B=False):
        if atom_A is False:
            try:
                atom_A = int(input('Please specify reference atom in first monomer: '))
                atom_B = int(input('Please specify reference atom in second monomer: '))
            except:
                print('Using defualt atom')
                atom_A = 3
                atom_B = 6
        #else:
        #    pass
        return (atom_A,atom_B)

    def resize(self, n_configs=False,configs=False,  first_n_configs=False, monomers=False, dis_lower=False, dis_upper=False, dis_new_lower=False, dis_new_upper=False):
        import numpy as np
        configs_new = list()
        configs = self.configs_check(configs)

        if dis_new_lower is False:#This is a weak argument, can be improved in many levels.
            try:
                #a = input('')
                n_configs = float(input('The number of new configuration you want is : '))
                dis_lower = float(input('The original distance_min (Angstrom) you want is : '))
                dis_upper = float(input('The original distance_max (Angstrom) you want is: '))

                dis_new_lower = float(input('The new distance_new_min (Angstrom) you want is: '))
                dis_new_upper = float(input('The new distance_new_max (Angstrom) you want is: '))
            except:
                print('Using default dis boundaries')
                dis_lower = float(2)
                dis_upper = float(5)
                dis_new_lower = float(6)
                dis_new_upper = float(9)

        monomers = self.monomers(configs)
        (atom_A,atom_B) = self.chosen_atom()

        if n_configs is not False:
            n_configs = int(n_configs)
        else:
            n_configs = len(configs)

        configs_count = 0
        configs_ok_count = 0
        rand_pool = len(configs) + 1
        while configs_ok_count < n_configs:
            configs_count += 1
            #print(np.random.randint(1,rand_pool))
            config = configs[np.random.randint(1,rand_pool)-1]
            #self.prt(config)
            dis = self.distance(config,atom_A,atom_B)
            if dis_lower <= dis <= dis_upper:
                configs_ok_count += 1
                dis_new = np.random.uniform(dis_new_lower,dis_new_upper)
                config_new = self.translate(atom_A,atom_B,dis_new,config)
                configs_new.append(config_new)
                #print('{:5.2f} {:5.2f} {:5.2f} {:5.2f}'.format(self.distance(config,atom_A,atom_B),self.distance(config_new,atom_A,atom_B),dis_new,self.distance(config_new,atom_A,atom_B)-self.distance(config,atom_A,atom_B)))
                print('New distance ({:5d}): {:5.2f}'.format(configs_ok_count,self.distance(config_new,atom_A,atom_B)))
        print('{:d}/{:d} configurations are returned as list.'.format(len(configs_new),configs_count))
        print('*Add point: Expand finished.')
        return configs_new #List of configs that have just been expanded.

    def vector(self,config,atom_A,atom_B=False):
        """For a given configuration and two atoms A and B, give back the vector BA.

        """
        import numpy as np

        if atom_B is False:
            atom_A = int(atom_A) - 1
            return np.array(config[2][atom_A][1])
        else:
            atom_A = int(atom_A) - 1
            atom_B = int(atom_B) - 1
            atom_A = np.array(config[2][atom_A][1])
            atom_B = np.array(config[2][atom_B][1])
            vector21 = atom_B-atom_A
            return vector21

    def molden(self,configs=False):
        """Give the configs in list form, and show them using molden"""

        configs = self.configs_check(configs)
        self.write('plot.temp',configs)
        address = self.alias('molden')
        self.cl(address +' plot.temp')

    def cl(self,command):
        """self.cl is to interact with terminal and catch the output.


        #ip::string, command line as string input
        #op::string, return value is the output of command line
        #Notice, each time when change dire.ctly, cl starts from currect directory.
        #Use three ' if you want to input multiple line
        """
        import subprocess
        import os
        import shlex
        arg = shlex.split(command)
        p = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
        (output, err) = p.communicate()
        #print(output)
        return output #output is bytes

    def alias(self,name,path=False):
        """Search ~/.bash_profile and return the alias cotent as a string.

        *Can overwrite default search path.
        """

        import re

        if path is False:
            path = '~/.bash_profile'
        address = 'Alias not found'
        #address_book = self.cl('cat '+path).decode("utf-8")
        address_book = self.cl('cat '+path)
        #print(type(address_book))
        line_count = 0
        for line in address_book.strip().split():
            line = line.decode('utf-8') #The line here could be a bytes type, so have to decode it.
            #print(line)
            try:
                if line[0:len(name)+1] == name+'=':  #different
                    expansion = line.split('"')[1]
                    break #If cannot find in single quote, then find double quote
            except:
                if line[0:len(name)+1] == name+'=':  # 'is' and '==' are different here
                    expansion = line.split("'")[1]
        #print(expansion)

        return expansion #This is expansion of alias

    def plot(self,configs = False,binwidth=False,ref=False):
        import numpy as np
        import matplotlib.pyplot as plt

        configs = self.configs_check(configs)

        if configs is not False:
            energy_array = list()
            for config in configs:
                energy_array.append(config[1][0][0]* self.hartree_to_cm)
            #print(energy_array)
        else:
            energy_array = self.energy_array_cm

        if ref is True:
            energy_array = energy_array-energy_array[0]
        fig = plt.figure()
        ax = fig.add_subplot(111)

        x = energy_array
        E_min = energy_array[0]
        E_max = energy_array[-1]
        #x = np.random.normal(0,1,1000)
        #print(x)

        if binwidth is False:
            binwidth = 50
        else:
            binwidth = int(binwidth)

        print('\nCurrent binwidth is {:d} cm-1\n'.format(binwidth))


        numBins = (self.energy_highest_cm - self.energy_lowest_cm) // binwidth
        #print(numBins)
        #numBins = 100
        ax.hist(x,numBins,color='green',alpha=0.8)
        #ax.boxplot(x,numBins)#,color='green',alpha=0.8)
        (count,x_tic) = np.histogram(x,numBins)


        count.sort()
        #print(count)
        count_highest = count[-1]
        ax.set_xlabel("Energy(cm$^{-1}$)")
        ax.set_ylabel("Count")
        #ax.set_title("Original data")
        ax.annotate('E_min:\n{:10.2f}'.format(E_min),xy=(E_min,0),xytext=(E_min, count_highest*0.2),arrowprops=dict(arrowstyle="->"))
        ax.annotate('E_max:\n{:10.2f}'.format(E_max),xy=(E_max,0),xytext=(E_max-(E_max-E_min)*0.1, count_highest*0.2),arrowprops=dict(arrowstyle="->"))
        #ax.xlabels[-1] = '300+'
        fig = plt.gcf()
        plt.show()
        decision = input("Do you want to save the file? (Enter 'y' to save, enter others to skip)")
        if decision is 'y':
            filename = input('Please specify .eps (1200 dpi) filename: ').strip()

            fig.savefig(filename, format='eps', dpi=1200)
            print('Plot saved to {}.'.format(filename))

    def plot2(self,clip_rate=False, binwidth=False,configs = False):
        import numpy as np
        import matplotlib.pyplot as plt

        configs = self.configs_check(configs)

        if clip_rate is False:
            clip_rate = 95
            print('Using default (95%) clip rate')
        clip_out = (100.0-float(clip_rate))/200
        print(clip_out)

        if configs is not False:
            energy_array = list()
            for config in configs:
                energy_array.append(config[1][0][0]* self.hartree_to_cm)
            #print(energy_array)
        else:
            energy_array = self.energy_array_cm
        fig = plt.figure()

        index = int(len(energy_array)*clip_out)

        print(index)
        energy_array = np.clip(energy_array,energy_array[index],energy_array[-index])
        E_min = energy_array[0]
        E_max = energy_array[-1]
        ax = fig.add_subplot(111)
        x = energy_array

        if binwidth is False:
            binwidth = 50
        else:
            binwidth = int(binwidth)

        print('\nCurrent binwidth is {:d} cm-1\n'.format(binwidth))

        fig.canvas.draw()

        numBins = (self.energy_highest_cm - self.energy_lowest_cm) // binwidth
        ax.hist(x,numBins,color='green',alpha=0.8)
        (count,x_tic) = np.histogram(x,numBins)
        count.sort()
        #print(count)
        #print(x_tic)
        count_highest = count[-1]
        ax.set_xlabel("Energy(cm$^{-1}$)")
        ax.set_ylabel("Count")
        ax.set_title("{:3.2f}% clipped distribution\nTotal count: {:d}\n".format(clip_rate,len(energy_array)))

        ax.annotate('\nTail {:3.2f}%\n< {:7.1f}cm$^{{-1}}$'.format(clip_out*100,E_min),xy=(E_min,0),xytext=(E_min-(E_max-E_min)*0.05, count_highest*0.2),arrowprops=dict(arrowstyle="->"))
        ax.annotate('\nHead {:3.2f}%\n> {:7.1f}cm$^{{-1}}$'.format(clip_out*100,E_max),xy=(E_max,0),xytext=(E_max+(E_max-E_min)*0.05, count_highest*0.2),arrowprops=dict(arrowstyle="->"))
        #ax.xlabels[-1] = '300+'
        #labels = ax.get_xticks().tolist()#To change specific labels
        #print(labels)
        # labels_count_min = 0
        # labels_count_max = 0
        #
        # for label in labels:
        #     print(label)
        #     if float(label) < E_min:
        #         labels[labels_count_min]=''
        #         labels_count_min += 1
        #
        #
        #     if float(label) < E_max:
        #         labels_count_max +=1
        #     else:
        #         labels[labels_count_max]=''
        #labels[labels_count_min + 1]=''
        #labels[labels_count_min] = '< {:<9.1f}'.format(float(int(E_min)))
        #labels[labels_count_max] = '> {:<9.1f}'.format(float(int(E_max)))
        #label()
        #labels[0]=''
        #labels[-1]=''
        #labels[1] = '< {:<9.1f}'.format(float(int(E_min)))
        #labels[-2]='> {:<9.1f}'.format(float(int(E_max)))
        #ax.set_xticklabels(labels)
        fig = plt.gcf()
        plt.show()
        try:
            decision = input("Do you want to save the file? (Enter 'y' to save, enter others to skip)")
            if decision is 'y':
                filename = input('Please specify .eps (1200 dpi) filename: ').strip()

                fig.savefig(filename, format='eps', dpi=1200)
                print('Plot saved to {}.'.format(filename))
        except:
            pass

    def v2b_configs(self,configs=False):
        """It returns configs of two monomers"""

        configs = self.configs_check(configs)
        monomers = self.monomers(configs)
        monomer_A = monomers[0]

        configs_A = list()
        molecule_count_A = len(monomer_A)

        for config in configs:
            coordinate=list()
            for atom in monomer_A:
                coordinate.append(config[2][atom-1])
            configs_A.append([[molecule_count_A],[[0]],coordinate])


        monomer_B = monomers[1]
        configs_B = list()
        molecule_count_B = len(monomer_B)
        for config in configs:
            coordinate=list()
            for atom in monomer_B:
                coordinate.append(config[2][atom-1])
            configs_B.append([[molecule_count_B],[[0]],coordinate])

        self.prt(configs_A)
        self.prt(configs_B)

    def int2list(self,str):
        """Convert the integer input to string into a list of numbers"""
        try:
            lst = list()
            for line in str.split():
                #print('TTTTest-----------------',line)
                if line is "\'":
                    continue
                lst.append(int(line.strip()))
        except:
            lst = str

        return lst #This list of int numbers

    def slice(self,monomer,configs=False):
        """Given the number of atoms, and split the molecule into list
            configs::given configuration
            monomers::the order of the monomer
        """
        configs = self.configs_check(configs)
        monomer = self.int2list(monomer)
        configs_monomer = list()
        for config in configs:
            config_monomer = list()
            config_monomer.append([len(monomer)])
            config_monomer.append(config[1])

            temp = list()
            for i in monomer:
                #print(config[2][i-1])
                temp.append(config[2][i-1])
            config_monomer.append(temp)
            #print(config_monomer)
            #self.prt(config_monomer)
            configs_monomer.append(config_monomer)
        #self.prt(configs_monomer[0])

        return configs_monomer

    def nodelist(self,configs=False):

        configs = self.configs_check(configs)
        cores = 24
        coresperjob = 4
        subjobs = 24
        nconfigs = len(configs)
        n = int(nconfigs / subjobs)
        lst = list()
        for i in range(1, subjobs + 1):
            start = (i - 1) * n
            end = start + n
            if i == subjobs:
                end = nconfigs
            lst.append((end - start))


        return lst

    def split2node(self,name,configs=False):

        configs=self.configs_check(configs)
        cores = 24
        coresperjob = 4
        subjobs = 24
        nconfigs = len(configs)
        n = int(nconfigs/subjobs)
        lst = self.nodelist(configs)
        for i in range(1,subjobs+1):
            start = (i-1)*n
            end = start+n
           # print('Start:{:d},end:{:d}'.format(start,end))
            if i == subjobs:
                end = nconfigs
               # print('This was run..................')
            #print('Start:{:d} End:{:d}'.format(start+1,end+1))
            #self.write('{}_{:02d}'.format(name,i),configs[start:end])
            #print('test2',configs[start])
            self.write(name,configs[start:end])
            self.cl("""mkdir {}{:02d}
                       mv {} {}{:02d}""".format(name,i,name,name,i))

        return lst #Describes the number of configs

    def pbs(self,ndata=1,usr='kee',job='spe'):
        """This is to take configurations and do ab initio calculation automatically"""


        #configs = self.configs_check(configs)
        #natom = configs[0][0]
        #ndata = len(configs)
        f = open('1_submit_sub.que','w')
        f.write('''#!/bin/bash
            #PBS -q xeon16
            #PBS -l nodes=1:ppn=4
            #PBS -N '''+job+'''
            #PBS -r n
            #PBS -c n
            #PBS -m n
            #PBS -e /dev/null
            #PBS -o /dev/null
            #PBS -S /bin/sh

            id=`echo $PBS_JOBID |cut -d. -f1`
            basename=`echo $infile |sed 's/\.[a-zA-Z]*$//' `
            log="$PBS_O_WORKDIR/spe.${id}.log"
            exe=${exe:="1_submit_sub.csh"}
            ndata="'''+str(ndata)+'''"

            export TMPDIR="/scratch/$USER"
            echo "job id:        " $PBS_JOBID  >> $log
            echo "input file:    " $infile $PBS_JOBNAME >> $log
            echo "exicutalbe:" $exe >> $log
            echo "job starts at: " `date` >> $log
            echo "submitted from:" $PBS_O_HOST >> $log
            echo "submitted to:  "  >> $log
            cat $PBS_NODEFILE >> $log
            echo "" >> $log

            MYDIR="/scratch/'''+usr+'''/spe-$id"
            mkdir $MYDIR
            cd $MYDIR
            cp $PBS_O_WORKDIR/* . -r

            echo ./$exe $ndata >>$log 2>&1
            ./$exe $ndata >>$log 2>&1
            cp -r $MYDIR/* $PBS_O_WORKDIR/
            rm -rf $MYDIR

            echo "job ends at:   " `date`  >> $log
                    ''')
        f.close()

        return None

    def molpro(self,file='b',natom=1,ndata=1,basis='avtz',ab='ccsd(t)-f12',title='co2h2o'):

        #configs = self.configs_check(configs)
        #natom = configs[0][0]
        #ndata = len(configs)
        #print(type(natom))
        #print('This is natom:',natom)
        #print(type(ab))
        #natom=6
        f = open('1_submit_sub.csh','w')
        #f.write('This is test')
        f.write('''#!/bin/csh -f
            #This file split a file into nfile files, and creat standard molpro input file for each configuration.


            #-----------------------User define-----------------
            set nfile = {0:d} #No. configs
            set inputfile = '{5}' #No suffix
            #---------------------------------------------------

            set n = {1:d}
            set i = 1
            while ( $i <= $nfile )
              set tail = $i
              if ( $i < "10" ) then
                set tail = 00$i
              else if ($i < "100") then
                set tail = 0$i
              endif

              set k = `expr $i - 0`
              set j = `expr $k \* $n`
              set input = $inputfile$tail


            echo '*** {2}'        >> $input.temp.head
            echo 'memory,100,m'      >> $input.temp.head
            echo 'basis {3} '        >> $input.temp.head
            echo 'geomtype = xyz'    >> $input.temp.head
            echo 'geometry = {{'      >> $input.temp.head


            #change here accroding to different input
            head -n $j $inputfile | tail -$n  >> $input.temp.config
            cat $input.temp.head  $input.temp.config >> $input
            rm $input.temp.head
            rm $input.temp.config

            echo ''                  >> $input
            echo '}}'                 >> $input
            echo 'hf'                >> $input
            echo '{4}'       >> $input
            echo '---'               >> $input

            /usr/local/bin/molprop_2010_1_Linux_x86_64_i8 -n 4 $input
              set  i = `expr $i + 1`
            end
        '''.format(ndata,natom+2,title,basis,ab,file))
        f.close()
        self.cl("chmod +x '1_submit_sub.csh' ")
        return None

    def submit(self,configs=False,file='dimer',v2b=False,basis='avtz',ab='ccsd(t)-f12',usr='kee',monomerA=None,monomerB=None):

        configs = self.configs_check(configs)
        nsubjobs=24
        molpro = '1_submit_sub.csh'
        pbs = '1_submit_sub.que'

        self.split2node('dimer',configs) #It splits file into nodes and return number of configs as list

        lst = self.nodelist(configs)
        natom = configs[0][0][0]
        #print(type(natom))



        if v2b is False:
            ans = input('Do you want two-body energy:? y/n \n')
            if ans is 'y':
                v2b = True


        if v2b is True:
            if monomerA is None:
                self.order(configs)
                monomerA = input('First monomer: ')
                self.order(configs)
                monomerB = input('Second monomer: ')

            monomerA = self.slice(self.int2list(monomerA),configs)
            monomerB = self.slice(self.int2list(monomerB),configs)


            self.split2node('monomerA',monomerA)
            self.split2node('monomerB',monomerB)
            #MonomerA
            natomA = monomerA[0][0][0]
            natomB = monomerB[0][0][0]

            #Generate molpro file and mv it to the dir

            for i in range(1,nsubjobs+1):
                ndata = lst[i-1]
                #print('------t1')

                self.molpro(file='monomerA',natom=natomA,ndata=ndata)
                self.cl('mv {} monomerA{:02}'.format(molpro, i))
                self.molpro(file='monomerB',natom=natomB,ndata=ndata)
                self.cl('mv {} monomerB{:02}'.format(molpro, i))
                self.pbs(ndata=ndata,job='monomerA')
                self.cl('mv {} monomerA{:02}'.format(pbs, i))
                self.pbs(ndata=ndata,job='monomerB')
                self.cl('mv {} monomerB{:02}'.format(pbs, i))
                self.cl('''cd monomerA{:02}
                            qsub {}'''.format(i,pbs))
                print(self.cl('''cd monomerB{:02}
                            qsub {} '''.format(i,pbs)))


        for i in range(1, nsubjobs + 1):
            ndata = lst[i - 1]
            self.molpro(file='dimer', natom=natom, ndata=ndata)
            self.cl('mv {} dimer{:02}'.format(molpro, i))
            self.pbs(ndata=ndata,job='dimer')
            self.cl('mv {} dimer{:02}'.format(pbs, i))
            self.cl('''cd dimer{:02}
                        qsub {} '''.format(i, pbs))



             #   self.cl('''mv dim
              #  ''')

        return None

    #def collect(self):
    #    nsubjobs = 24

     #   for i in range(1, nsubjobs + 1):

    def extract_ee(self,info=None):

        f = open(dimer001.out)
        i = 0
        lst = list()
        for line in f:
            line = line.strip()
            line = line.split()
            i += 1
            if len(line) ==0:
                continue
            if line[0] is 'geometry':
                lst.append(i)  #Record the line
        return

    def extract(self,configs=False,v2b=False,monomerA=None,monomerB=None):

        configs = self.configs_check(configs)
        nsubjobs = 24
        molpro = '1_submit_sub.csh'
        pbs = '1_submit_sub.que'
        extract = '2_extractE_sub.csh'
        print('Warning!!!!Use unsorted list to do extract')

        lst = self.nodelist(configs)  # It splits file into nodes and return number of configs as list
        print(lst)

        natom = configs[0][0][0]
        # print(type(natom))

        if v2b is True:
            self.cl('rm monomerA.abE')
            self.cl('rm monomerB.abE')
            if monomerA is None:
                self.order(configs)
                monomerA = input('First monomer: ')
                self.order(configs)
                monomerB = input('Second monomer: ')

            monomerA = self.int2list(monomerA)
            monomerB = self.int2list(monomerB)


            # MonomerA
            natomA = len(monomerA)
            natomB = len(monomerB)

            # Generate molpro file and mv it to the dir

            for i in range(1, nsubjobs + 1):
                ndata = lst[i - 1]
                print('Extracting monomerA {:d}'.format(i))
                #print('natomA: {:d}'.format(natomA))
                #print(monomerA)
                #print(len(monomerA))
                #for el in monomerA:
                #    print (el)
                self.extract_monomer(file='monomerA',nfile=ndata,natom=natomA)
                self.cl('mv {} monomerA{:02}'.format(extract, i))
                print('Extracting monomerB {:d}'.format(i))
                self.extract_monomer(file='monomerB', nfile=ndata, natom=natomB)

                self.cl('mv {} monomerB{:02}'.format(extract, i))

                self.cl('''cd monomerA{:02}
                            rm monomerA.abE
                                     ./{}
                            cp monomerA.abE monomerA.abE{:02d}
                            mv monomerA.abE{:02d} ../
                            cd ../
                            cat monomerA.abE{:02d} >> monomerA.abE
                            rm monomerA.abE{:02d}'''.format(i, extract,i,i,i,i))
                self.cl('''cd monomerB{:02}
                            rm monomerB.abE
                                     ./{}
                            cp monomerB.abE monomerB.abE{:02d}
                            mv monomerB.abE{:02d} ../
                            cd ../
                            cat monomerB.abE{:02d} >> monomerB.abE
                            rm monomerB.abE{:02d}'''.format(i, extract, i, i, i, i))


        self.cl('rm dimer.abE')
        for i in range(1, nsubjobs + 1):
            print('Extracting dimer {:d}'.format(i))
            ndata = lst[i - 1]
            self.extract_dimer(file='dimer', nfile=ndata, natom=natom)
            self.cl('mv {} dimer{:02}'.format(extract, i))
            self.cl('''cd dimer{:02}
                        rm dimer.abE
                                 ./{}
                        cp dimer.abE dimer.abE{:02d}
                        mv dimer.abE{:02d} ../
                        cd ../
                        cat dimer.abE{:02d} >> dimer.abE
                        rm dimer.abE{:02d}'''.format(i, extract, i, i, i, i))

        return None

    def extract_dimer(self,file='dimer',nfile=15,natom=6):


        f2 = open('2_extractE_sub.csh','w')
        f2.write('''#!/bin/csh -f
            #When using this file, should be used for molpro outputfiles of the same batch. namexxx.out, and should be in same output format. Only the fist colum of $output will be extracted.
            #Need to set the number of output files need to be processed.
            #Need to set the number of atoms of the batch config
            #Need to set the order of atoms, i.e. where is the line of atom coordinates.
            #nfile=number of output files in this folder

            #please change following settings accordingly.

            #----------------------------User define-----------------
            set input = '{}'
            set nfile = {:d}   #Num of file per node
            #--------------------------------------------------------


            set output = '{}.abE'
            set natom = {:d}
            set lineAtom1 = 29
            set lineAtom2 = 30
            set lineAtom3 = 31
            set lineAtom4 = 32
            set lineAtom5 = 33
            set lineAtom6 = 34

            set i = 1
            while ( $i <= $nfile )
              set tail = $i
              if ( $i < "10" ) then
                set tail = 00$i
              else if ($i < "100") then
                set tail = 0$i
              endif

            #Make sure the line number of $output and coordinate are corret, using keyword1 and keyword2.
            #if ((head -n $keyword1 $input$tail.out | tail -1) == $linekey1 && (head -n $keyword2 $input$tail.out | tail -1) == $linekey2)


            #Number of atom
            echo $natom >> $output

            #Energy extracted from output files.

            tail -n 3 $input$tail.out |  head -1 | awk '{{print $1}}' >> $output



            ####What arrangement of atoms do you want?
            head -n $lineAtom1 $input$tail.out | tail -1 >> $output
            head -n $lineAtom2 $input$tail.out | tail -1 >> $output
            head -n $lineAtom3 $input$tail.out | tail -1 >> $output
            head -n $lineAtom4 $input$tail.out | tail -1 >> $output
            head -n $lineAtom5 $input$tail.out | tail -1 >> $output
            head -n $lineAtom6 $input$tail.out | tail -1 >> $output
            #head -n $lineAtom7 $input$tail.out | tail -1 >> $output
            #head -n $lineAtom8 $input$tail.out | tail -1 >> $output
            #head -n $lineAtom9 $input$tail.out | tail -1 >> $output
            #head -n $lineAtom10 $input$tail.out | tail -1 >> $output

            set i = `expr $i + 1`

            #else
            #echo 'There is a bad file in file $i' >> badfilereport
            #endif
            end

             '''.format(file,nfile,file,natom))
        f2.close()
        self.cl('chmod +x 2_extractE_sub.csh')

        return None

    def extract_monomer(self, file='monomer', nfile=15, natom=6):

        f2 = open('2_extractE_sub.csh', 'w')
        f2.write('''#!/bin/csh -f
               #When using this file, should be used for molpro outputfiles of the same batch. namexxx.out, and should be in same output format. Only the fist colum of $output will be extracted.
               #Need to set the number of output files need to be processed.
               #Need to set the number of atoms of the batch config
               #Need to set the order of atoms, i.e. where is the line of atom coordinates.
               #nfile=number of output files in this folder

               #please change following settings accordingly.

               #----------------------------User define-----------------
               set input = '{}'
               set nfile = {:d}   #Num of file per node
               #--------------------------------------------------------


               set output = '{}.abE'
               set natom = {:d}
               set lineAtom1 = 29
               set lineAtom2 = 30
               set lineAtom3 = 31
               set lineAtom4 = 32
               set lineAtom5 = 33
               set lineAtom6 = 34

               set i = 1
               while ( $i <= $nfile )
                 set tail = $i
                 if ( $i < "10" ) then
                   set tail = 00$i
                 else if ($i < "100") then
                   set tail = 0$i
                 endif

               #Make sure the line number of $output and coordinate are corret, using keyword1 and keyword2.
               #if ((head -n $keyword1 $input$tail.out | tail -1) == $linekey1 && (head -n $keyword2 $input$tail.out | tail -1) == $linekey2)


               #Number of atom
               echo $natom >> $output

               #Energy extracted from output files.

               tail -n 3 $input$tail.out |  head -1 | awk '{{print $1}}' >> $output



               ####What arrangement of atoms do you want?
               head -n $lineAtom1 $input$tail.out | tail -1 >> $output
               head -n $lineAtom2 $input$tail.out | tail -1 >> $output
               head -n $lineAtom3 $input$tail.out | tail -1 >> $output
               #head -n $lineAtom4 $input$tail.out | tail -1 >> $output
               #head -n $lineAtom5 $input$tail.out | tail -1 >> $output
               #head -n $lineAtom6 $input$tail.out | tail -1 >> $output
               #head -n $lineAtom7 $input$tail.out | tail -1 >> $output
               #head -n $lineAtom8 $input$tail.out | tail -1 >> $output
               #head -n $lineAtom9 $input$tail.out | tail -1 >> $output
               #head -n $lineAtom10 $input$tail.out | tail -1 >> $output

               set i = `expr $i + 1`

               #else
               #echo 'There is a bad file in file $i' >> badfilereport
               #endif
               end

                '''.format(file, nfile, file, natom))
        f2.close()
        self.cl('chmod +x 2_extractE_sub.csh')

        return None

    def v2b(self,configs=False,configs1=False,configs2=False):
        configs = self.configs_check(configs)
        configs1 = self.configs_check(configs1)
        configs2 = self.configs_check(configs2)


        i = 0
        for config in configs:
            #print(config[1][0][0],'-',configs1[i][1][0][0],'-',configs2[i][1][0][0])
            config[1][0][0] = config[1][0][0]- configs1[i][1][0][0] - configs2[i][1][0][0]
            #print(config[1][0][0])
            i += 1

        return configs

    def compare(self,configs1=False,configs2=False,configs3=False,atomA=1,atomB=2):

        import numpy as np
        import matplotlib.pyplot as plt

        aucm = 219474.63
        configs1 = self.configs_check(configs1)
        #configs1 = self.sort(configs1)
        configs2 = self.configs_check(configs2)
        configs3 = self.configs_check(configs3)
        dis1 = list()
        dis2 = list()
        dis3 = list()
        e1 = list()
        e2 = list()
        e3 = list()
        ecompare = list()
        i = 0
        for config in configs1:

            dis1.append(self.distance(config,atomA,atomB))
            dis2.append(self.distance(configs2[i],atomA,atomB))
            dis3.append(self.distance(configs3[i],atomA,atomB))
            #print(dis2[i])
            #self.prt(configs2[i])
            e1.append(config[1][0][0]*aucm)
            e2.append(configs2[i][1][0][0]*aucm)
            e3.append(configs3[i][1][0][0] * aucm)
            ecompare.append(e3[i]-e2[i])
            i=i+1
        s = 7
        fig = plt.figure()
        ax = fig.add_subplot(211)
        ax.set_title('Long-range-fit')
        ax.scatter(dis1,e1,s=s,color='g',label='a0=4')
        #ax.annotate('E_min:\n{:10.2f}'.format(E_min), xy=(E_min, 0), xytext=(E_min, count_highest * 0.2),
                    #arrowprops=dict(arrowstyle="->"))
        ax.scatter(dis2,e2,s=s,color='r',label='a0=5')
        ax.scatter(dis3,e3,s=s,color='b',label='a0=6')
        #ax.scatter(dis3,e3,s=5, color='b', label='a0=6')
        ax.legend()
        ax2 = fig.add_subplot(212)
        ax2.scatter(dis1, ecompare, color='r', label='(whole range) - (ab initio)')
        ax2.legend()
        print('show')
        plt.show()


        return None

    def compare_cm(self,configs1=False,atomA=3,atomB=6):

        import numpy as np
        import matplotlib.pyplot as plt

        aucm = 219474.63

        configs1 = self.configs_check(configs1)
        configs1 = self.sort(configs1)

        dis1 = list()

        e1 = list()

        i = 0
        for config in configs1:

            dis1.append(self.distance(config,atomA,atomB))
            e1.append(config[1][0][0]*aucm)


        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.scatter(dis1,e1)
        plt.show()


        return None

    def compared(self,configss, atomA=3,atomB=6,s=100, title='',color='b',marker='.',label='default',xmin=None,xmax=None,ymin=None,ymax=None,xmin2=None,xmax2=None,ymin2=None,ymax2=None,loc=None,fontsize=12):

        import numpy as np
        import matplotlib.pyplot as plt

        aucm = 219474.63

        n = len(configss)
        dis = list()
        e = list()
        ecompare=list()
        for i in range(0,n):
            configss[i] = self.configs_check(configss[i])
            configss[i] = self.sort(configss[i], key='distance', subkey1=atomA, subkey2=atomB)



        fig = plt.figure(figsize=(10,6 ))
        ax = fig.add_subplot(211)
        axes = plt.gca()
        axes.set_xlim([xmin, xmax])
        plt.ylim(ymin, ymax)
        ax2 = fig.add_subplot(212)#,sharex=True)
        axes = plt.gca()
        axes.set_xlim([xmin2, xmax2])
        plt.ylim(ymin2, ymax2)

        ax.set_title(title)


        for n in range(0,n):
            print(n)
            dis_temp = list()
            e_temp = list()
            e_temp_ref = list()
            i = 0
            j = 0
            for config in configss[n]:


                if j%1 == 0 :
                    #j = j + 1
                    print(j)
                    print('OK')
                    print('i=',i)
                    #continue

                    dis_temp.append(self.distance(config,atomA,atomB))
                    e_temp.append(config[1][0][0]*aucm)
                    e_temp_ref.append(configss[-1][j][1][0][0]*aucm-e_temp[i])
                    i=i+1
                j=j+1
            ecompare.append(e_temp_ref)
            dis.append(dis_temp)
            e.append(e_temp)
            #print(dis[n])
            #ax.plot(dis[n], e[n], s=s, color=color[n],marker=marker[n], label=label[n],facecolors='none',linewidth='2')

            ax.plot(dis[n], e[n], ) #color=color[n],marker=marker[n], label=label[n],facecolors='none',linewidth='2')
            ax2.scatter(dis[n], ecompare[n], s=s, color=color[n],marker=marker[n], label=label[n],facecolors='none',linewidth='2')


        #ax.scatter(dis1,e1,s=s,color='g',label='a0=4')
        #ax.annotate('E_min:\n{:10.2f}'.format(E_min), xy=(E_min, 0), xytext=(E_min, count_highest * 0.2),
                    #arrowprops=dict(arrowstyle="->"))
        #ax.scatter(dis2,e2,s=s,color='r',label='a0=5')
        #ax.scatter(dis3,e3,s=s,color='b',label='a0=6')
        #ax.scatter(dis3,e3,s=5, color='b', label='a0=6')
        #plt.ylim(-10, 0.1)

        ax.legend(loc=4,fontsize=12)
        ax2.legend(loc=4, fontsize=12)
        #ax.set_xlabel('$r_5$ ($\AA$)')
        ax.set_ylabel('V$_{\mathrm{2b}}$ (cm$^{-1}$)',fontsize=fontsize)
        ax.text(4.1,-5,'(a)',fontsize=fontsize)
        ax2.text(4.1, 5.8, '(b)',fontsize=fontsize)
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(fontsize)

        for tick in ax2.xaxis.get_major_ticks():
            tick.label.set_fontsize(fontsize)


        ax2.set_xlabel('$r_5$ ($\AA$)',fontsize=fontsize)
        ax2.set_ylabel('$\Delta \mathrm{V}_{\mathrm{2b}}$ (cm$^{-1}$)',fontsize=fontsize)


        plt.show()
        #decision = input("Do you want to save the file? (Enter 'y' to save, enter others to skip)")
        decision = 'y'
        if decision is 'y':
            #filename = input('Please specify .eps (1200 dpi) filename: ').strip()
            filename = 'switch.eps'
            fig.savefig(filename, format='eps', dpi=1200)
            print('Plot saved to {}.'.format(filename))



        return None

    def plotv2b(self, configs=False, binwidth=False):
        import numpy as np
        import matplotlib.pyplot as plt



        SMALL_SIZE = 12
        MEDIUM_SIZE = 14
        BIGGER_SIZE = 16

        plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
        plt.rc('axes', titlesize=SMALL_SIZE)  # fontsize of the axes title
        plt.rc('axes', labelsize=BIGGER_SIZE)  # fontsize of the x and y labels
        plt.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
        plt.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
        plt.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
        plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


        configs = self.configs_check(configs)

        if configs is not False:
            energy_array = list()
            for config in configs:
                energy_array.append(config[1][0][0] * self.hartree_to_cm)
                # print(energy_array)
        else:
            energy_array = self.energy_array_cm
        fig = plt.figure(figsize=(6,4))
        #plt.gcf().subplots_adjust(bottom=0.15)
        ax = fig.add_subplot(111)
        axes = plt.gca()
        axes.set_xlim([-1500, 3100])
        #plt.ylim(0, ymax)

        x = energy_array
        E_min = energy_array[0]
        E_max = energy_array[-1]
        # x = np.random.normal(0,1,1000)
        # print(x)

        if binwidth is False:
            binwidth = 50
        else:
            binwidth = int(binwidth)

        print('\nCurrent binwidth is {:d} cm-1\n'.format(binwidth))

        numBins = (self.energy_highest_cm - self.energy_lowest_cm) // binwidth
        # print(numBins)
        # numBins = 100
        ax.hist(x, numBins, color='white', alpha=0.5)
        # ax.boxplot(x,numBins)#,color='green',alpha=0.8)
        (count, x_tic) = np.histogram(x, numBins)

        count.sort()

        count_highest = count[-1]
        ax.set_xlabel("V$_{\mathrm{2b}}$ (cm$^{-1}$)")
        ax.set_ylabel("Count")

        fig = plt.gcf()
        plt.tight_layout()
        plt.show()
        #decision = input("Do you want to save the file? (Enter 'y' to save, enter others to skip)")
        decision = 'y'
        if decision is 'y':
            filename = 'v2b-energy-distribution.eps'#input('Please specify .eps (1200 dpi) filename: ').strip()

            fig.savefig(filename, format='eps', dpi=1200)
            print('Plot saved to {}.'.format(filename))

    def plotevsr_for_publication(self, configs1,configs2, atom_A=3, atom_B = 6, binwidth=False):
        "Plot the energy versus r_co distribution"
        import numpy as np
        import matplotlib.pyplot as plt



        SMALL_SIZE = 12
        MEDIUM_SIZE = 14
        BIGGER_SIZE = 16

        plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
        plt.rc('axes', titlesize=SMALL_SIZE)  # fontsize of the axes title
        plt.rc('axes', labelsize=BIGGER_SIZE)  # fontsize of the x and y labels
        plt.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
        plt.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
        plt.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
        plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


        configs1 = self.configs_check(configs1)
        configs2 = self.configs_check(configs2)


        energy_array1 = list()
        energy_array2 = list()
        distance_array1 = list()
        distance_array2 = list()
        for config in configs1:
            energy_array1.append(config[1][0][0] * self.hartree_to_cm)
            distance_array1.append(self.distance(config,atom_A,atom_B))
        for config in configs2:
            energy_array2.append(config[1][0][0] * self.hartree_to_cm)
            distance_array2.append(self.distance(config,atom_A,atom_B))

        fig = plt.figure(figsize=(6,4))
        #plt.gcf().subplots_adjust(bottom=0.15)
        ax1 = fig.add_subplot(111)
        #ax2 = fig.add_subplot(212)
        axes = plt.gca()
        axes.set_xlim([1.5, 17])
        axes.set_ylim([-2000, 4300])
        #axes.set_xlim([5, 24])
        #axes.set_ylim([-100, 50])
        #plt.ylim(0, ymax)


        # x = np.random.normal(0,1,1000)
        # print(x)
        ax1.scatter(distance_array1,energy_array1, s=1, c='k',marker='.',linewidths=0,edgecolors=None)
        left, bottom, width, height = [0.42, 0.52, 0.50, 0.4]
        ax2 = fig.add_axes([left, bottom, width, height])
        axes = plt.gca()
        axes.set_xlim([5, 23])
        axes.set_ylim([-150, 60])
        ax2.scatter(distance_array2, energy_array2, s=1, c='m', marker='.', linewidths=0, edgecolors=None)


        ax1.set_xlabel(r'$R_{\mathrm{CO}}$ ($\mathrm{\AA}$)')
        ax1.set_ylabel(r'V$_{\mathrm{2b}}$  (cm$^{-1}$)')
        ax1.text(11, -1500, 'fit-SR database', fontsize=14)
        ax2.text(13, -120, 'fit-LR database', fontsize=14,color='m')
        fig = plt.gcf()
        plt.tight_layout()
        plt.show()
        #decision = input("Do you want to save the file? (Enter 'y' to save, enter others to skip)")
        decision = 'y'
        if decision is 'y':
            filename = 'v2b-energy-vs-rco.eps'#input('Please specify .eps (1200 dpi) filename: ').strip()

            fig.savefig(filename, format='eps', dpi=1200)
            print('Plot saved to {}.'.format(filename))

    def pot_sum(self,configs=None):
        """Add all potentials (for v2b and v3b) for configs in clathrate"""
        configs = self.configs_check(configs)
        pot = 0
        for config in configs:
            pot = pot + config[1][0][0]
        print('Overall potential is {:14.8f}'.format(pot))
        return pot
    def pot_avg(self,configs=None):
        """Add all potentials (for v2b and v3b) for configs in clathrate"""
        configs = self.configs_check(configs)
        pot = 0
        for config in configs:
            pot = pot + config[1][0][0]
        print('Average potential is {:14.8f}'.format(pot/len(configs)))
        return pot/len(configs)
    def sort_unique(self, configs = False):
        """Sorted out all configs with unique potential"""
        configs = self.configs_check(configs)
        sort_configs = self.sort()
        uconfig = list()
        pot = 0

        for config in  sort_configs:
            if config[1][0][0] == pot:
                continue
            else:
                uconfig.append(config)
                pot = config[1][0][0]
        print('There are {:d} configs with unique values in all {:d} configs.'.format(len(uconfig),len(configs)))

        return uconfig




"""To keep this script as clean as possible, please use another script for test arguments"""

