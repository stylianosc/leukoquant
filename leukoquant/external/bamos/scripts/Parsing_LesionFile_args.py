
import os
import sys
import getopt

import re
import pandas as pd


def parse(filepath, numb_modal):
    """
    Parse text at given filepath

    Parameters
    ----------
    filepath : str
        Filepath for file to be parsed

    Returns
    -------
    data : pd.DataFrame
        Parsed data

    """

    data = []
    texture = []
    texture1 = []
    edginess = []
    mahalquant = []

    with open(filepath, 'r') as file:
        line = file.readline()
        while line:
            reg_match = _RegExLib(line)
            if reg_match.label:
                label = reg_match.label.group(1)
                print(label)
            if reg_match.ot:
                
                ot = reg_match.ot.group(1)
            if reg_match.ls:
                
                ls = reg_match.ls.group(1)
            if reg_match.mean_code:
                
                mean_code = reg_match.mean_code.group(1)
            if reg_match.volume:
                
                volume = reg_match.volume.group(1)
            if reg_match.surface:
                
                surface = reg_match.surface.group(1)
            if reg_match.sav:
                
                sav = reg_match.sav.group(1)
            if reg_match.edginess:
                edginess = []
                for i in range(0, numb_modal+1):
                    line = next(file)
                    value_list = line.strip().split(' ')
                    edginess.append(value_list[:-1])
            if reg_match.compactess:
                
                compactess = reg_match.compactess.group(1)
            if reg_match.mahalquant:
                mahalquant = []
                for i in range(0, numb_modal+1):
                    line = next(file)
                    value_list = line.strip().split(' ')
                    mahalquant.append(value_list[:-1])
            if reg_match.mahalmin:
                mahalmin_temp = reg_match.mahalmin.group(1)
                mahalmin = mahalmin_temp.strip().split(' ')
            if reg_match.mahalmean:
                mahalmean_temp = reg_match.mahalmean.group(1)
                mahalmean = mahalmean_temp.split(' ')
            if reg_match.mahalmax:
                mahalmax_temp = reg_match.mahalmax.group(1)
                mahalmax = mahalmax_temp.strip().split(' ')
            if reg_match.mean:
                mean_temp = reg_match.mean.group(1)
                mean = mean_temp.strip().split(' ')
            if reg_match.variance:
                variance_temp = reg_match.variance.group(1)
                variance = variance_temp.split(' ')
            if reg_match.distancegrav:
                distancegrav = reg_match.distancegrav.group(1)
            if reg_match.propneightwm:
                propneightwm = reg_match.propneightwm.group(1)
            if reg_match.propneighgm:
                propneighgm = reg_match.propneighgm.group(1)
            if reg_match.propneighcsf:
                propneighcsf = reg_match.propneighcsf.group(1)
            if reg_match.propneighnb:
                propneighnb = reg_match.propneighnb.group(1)
            if reg_match.propneighout:
                propneighout = reg_match.propneighout.group(1)
            if reg_match.grav:
                grav_temp = reg_match.grav.group(1)
                grav = grav_temp.split(' ')
            if reg_match.diffgrav:
                diffgrav_temp = reg_match.diffgrav.group(1)
                diffgrav = diffgrav_temp.split(' ')
            if reg_match.dgm:
                dgm = reg_match.dgm.group(1)
            if reg_match.sppot:
                sppot = reg_match.sppot.group(1)
            if reg_match.propicsf:
                propicsf = reg_match.propicsf.group(1)
            if reg_match.propecsf:
                propecsf = reg_match.propecsf.group(1)
            if reg_match.propcgm:
                propcgm = reg_match.propcgm.group(1)
            if reg_match.propdgm:
                propdgm = reg_match.propdgm.group(1)
            if reg_match.propwmh:
                propwmh = reg_match.propwmh.group(1)
            if reg_match.propartefact:
                propartefact = reg_match.propartefact.group(1)
            if reg_match.propneighartefact:
                propneighartefact = reg_match.propneighartefact.group(1)
            if reg_match.propborderextent:
                propborderextent_temp = reg_match.propborderextent.group(1)
                _ = propborderextent_temp.split(' ')  # parsed but not emitted (commented out in dict)
            if reg_match.distventr:
                distventr_temp = reg_match.distventr.group(1)
                distventr = distventr_temp.split(' ')
            if reg_match.distgmi:
                distgmi_temp = reg_match.distgmi.group(1)
                distgmi = distgmi_temp.split(' ')
            if reg_match.distwmi:
                distwmi_temp = reg_match.distwmi.group(1)
                distwmi = distwmi_temp.split(' ')
            if reg_match.distcsfi:
                distcsfi_temp = reg_match.distcsfi.group(1)
                distcsfi = distcsfi_temp.split(' ')
            if reg_match.distouti:
                distouti_temp = reg_match.distouti.group(1)
                distouti = distouti_temp.split(' ')
            if reg_match.distsp:
                distsp_temp = reg_match.distsp.group(1)
                distsp = distsp_temp.split(' ')
            if reg_match.ratiobb:
                ratiobb_temp = reg_match.ratiobb.group(1)
                ratiobb = ratiobb_temp.split(' ')
            if reg_match.bbrelgrav:
                bbrelgrav_temp = reg_match.bbrelgrav.group(1)
                bbrelgrav = bbrelgrav_temp.split(' ')
            if reg_match.extentbox:
                extentbox_temp = reg_match.extentbox.group(1)
                extentbox = extentbox_temp.split(' ')
            if reg_match.ratioext:
                ratioext_temp = reg_match.ratioext.group(1)
                ratioext = ratioext_temp.split(' ')
            if reg_match.eigenvalues:
                eigenvalues_temp = reg_match.eigenvalues.group(1)
                eigenvalues = eigenvalues_temp.split(' ')
            if reg_match.ratioeigen:
                ratioeigen_temp = reg_match.ratioeigen.group(1)
                ratioeigen = ratioeigen_temp.split(' ')
            if reg_match.proporigin:
                proporigin_temp = reg_match.proporigin.group(1)
                proporigin = proporigin_temp.split(' ')
            if reg_match.texture:
                texture.append(reg_match.texture.group(1))
            if reg_match.texture1:
                texture1.append(reg_match.texture1.group(1))

            if len(texture1) == numb_modal:
                    dict_of_data = {
                        'file': filepath,
                        'label': label,
                        'ot': ot,
                        'ls': ls,
                        'mean_code': mean_code,
                        'volume': volume,
                        'surface': surface,
                        'sav': sav,
                        'edginess': edginess,
                        'compactess': compactess,
                        'mahalquant': mahalquant,
                        # 'mahalmin': (mahalmin, ['MahalMin_%d' % i for i in
                        #                         range(0,numb_modal+1) ]),
                        # 'mahalmean': mahalmean,
                        # 'mahalmax': mahalmax,
                        # 'mean': mean,
                        # 'variance': variance,
                        'distancegrav': distancegrav,
                        'propneightwm': propneightwm,
                        'propneighgm': propneighgm,
                        'propneighcsf': propneighcsf,
                        'propneighnb': propneighnb,
                        'propneighout': propneighout,
                        'grav': grav,
                        # 'diffgrav': diffgrav,
                        'dgm': dgm,
                        'sppot': sppot,
                        'propicsf': propicsf,
                        'propecsf': propecsf,
                        'propcgm': propcgm,
                        'propdgm': propdgm,
                        'propwmh': propwmh,
                        'propartefact': propartefact,
                        'propneighartefact': propneighartefact,
                        # 'propborderextent': propborderextent,
                        # 'distventr': distventr,
                        # 'distgmi': distgmi,
                        # 'distwmi': distwmi,
                        # 'distcsfi': distcsfi,
                        # 'distouti': distouti,
                        # 'distsp': distsp,
                        # 'ratiobb': ratiobb,
                        # 'bbrelgrav': bbrelgrav,
                        # 'extentbox': extentbox,
                        # 'ratioext': ratioext,
                        # 'eigenvalues': eigenvalues,
                        # 'ratioeigen': ratioeigen,
                        # 'proporigin': proporigin,
                        'texture': texture,
                        'texture1': texture1,
                    }
                    for i in range(0, numb_modal+1):
                        dict_of_data['mahalmin%d' %i] = mahalmin[i]
                        dict_of_data['mahalmean%d' % i] = mahalmean[i]
                        dict_of_data['mahalmax%d' % i] = mahalmax[i]
                    for d in range(0, 3):
                        dict_of_data['bbrelgravmin%d' %d] = bbrelgrav[2*d]
                        dict_of_data['bbrelgravmax%d' % d] = bbrelgrav[2 *
                                                                          d+1]
                        dict_of_data['diffgrav%d' % d] = diffgrav[d]
                        dict_of_data['eigenvalues%d' % d] = eigenvalues[d]
                        dict_of_data['extentbox%d' % d] = extentbox[d]
                        dict_of_data['grav%d' % d] = grav[d]
                        dict_of_data['ratioeigen%d' % d] = ratioeigen[d]
                        dict_of_data['ratioext%d' % d] = ratioext[d]
                        dict_of_data['distcsfi%d' % d] = distcsfi[d]
                        dict_of_data['distgmi%d' % d] = distgmi[d]
                        dict_of_data['distwmi%d' % d] = distwmi[d]
                        dict_of_data['distouti%d' % d] = distouti[d]
                        dict_of_data['distventr%d' % d] = distventr[d]
                        dict_of_data['distsp%d' % d] = distsp[d]
                    for d in range(0, numb_modal):
                        dict_of_data['mean%d' % d] = mean[d]
                        dict_of_data['variance%d' % d] = variance[d]
                    for d in range(0, 4):
                        dict_of_data['ratiobb%d' % d] = ratiobb[d]
                        dict_of_data['proporigin%d' % d] = proporigin[d]



                    data.append(dict_of_data)
                    mahalquant=[]
                    edginess=[]
                    texture=[]
                    texture1=[]
            try:
                line = next(file)
            except Exception:
                break


    return data


class _RegExLib:
    """Set up regular expressions"""
    # use https://regexper.com to visualise these if required
    _reg_label = re.compile('Label (.*)\n')
    _reg_ot = re.compile('OutlierType (.*)\n')
    _reg_ls = re.compile('LesionSub (.*)\n')
    _reg_mean_code = re.compile('MeanCode (.*)\n')
    _reg_volume = re.compile('Volume (.*)\n')
    _reg_surface = re.compile('Surface (.*)\n')
    _reg_sav = re.compile('SAV (.*)\n')
    _reg_edginess = re.compile('Edginess (.*)\n')
    _reg_compactess = re.compile('Compactness (.*)\n')
    _reg_mahalquant = re.compile('MahalVecQuant (.*)\n')
    _reg_mahalmin = re.compile('MahalVecMin (.*)\n')
    _reg_mahalmean = re.compile('MahalVecMean (.*)\n')
    _reg_mahalmax = re.compile('MahalVecMax (.*)\n')
    _reg_mean = re.compile('MeanLesion (.*)\n')
    _reg_variance = re.compile('VarianceLesion (.*)\n')
    _reg_distancegrav = re.compile('DistanceToCentreGrav (.*)\n')
    _reg_propneightwm = re.compile('ProportionNeighWM (.*)\n')
    _reg_propneighgm = re.compile('ProportionNeighGM (.*)\n')
    _reg_propneighcsf = re.compile('ProportionNeighCSF (.*)\n')
    _reg_propneighnb = re.compile('ProportionNeighNB (.*)\n')
    _reg_propneighout = re.compile('ProportionNeighOut (.*)\n')
    _reg_grav = re.compile('CentreGravity (.*)\n')
    _reg_diffgrav = re.compile('DiffCentreGravity (.*)\n')
    _reg_dgm = re.compile('DGMBelonging (.*)\n')
    _reg_sppot = re.compile('SPPotProp (.*)\n')
    _reg_propicsf = re.compile('PropICSF (.*)\n')
    _reg_propecsf = re.compile('PropECSF (.*)\n')
    _reg_propcgm = re.compile('PropCGM (.*)\n')
    _reg_propdgm = re.compile('PropDGM (.*)\n')
    _reg_propwmh = re.compile('PropWMH (.*)\n')
    _reg_propartefact = re.compile('PropArtefact (.*)\n')
    _reg_propneighartefact = re.compile('PropNeighbourArtefact (.*)\n')
    _reg_propborderextent = re.compile('PropBorderExtent (.*)\n')
    _reg_distventr = re.compile('DistanceVentricle (.*)\n')
    _reg_distgmi = re.compile('DistanceGMI (.*)\n')
    _reg_distwmi = re.compile('DistanceWMI (.*)\n')
    _reg_distcsfi = re.compile('DistanceCSFI (.*)\n')
    _reg_distouti = re.compile('DistanceOutI (.*)\n')
    _reg_distsp = re.compile('DistanceSP (.*)\n')
    _reg_ratiobb = re.compile('RatioBB (.*)\n')
    _reg_bbrelgrav = re.compile('BoundingBoxRelGrav (.*)\n')
    _reg_extentbox = re.compile('ExtentBox (.*)\n')
    _reg_ratioext = re.compile('RatioExtentVolume (.*)\n')
    _reg_eigenvalues = re.compile('EigenValues (.*)\n')
    _reg_ratioeigen = re.compile('RatioEigen (.*)\n')
    _reg_proporigin = re.compile('ProportionOrigin (.*)\n')
    _reg_texture = re.compile('TextureDescriptors (.*)\n')
    _reg_texture1 = re.compile('TextureDescriptors (.*)\n')



    def __init__(self, line):
        # check whether line has a positive match with all of the regular expressions
        self.label = self._reg_label.match(line)
        self.ot = self._reg_ot.match(line)
        self.ls = self._reg_ls.match(line)
        self.mean_code = self._reg_mean_code.match(line)
        self.volume = self._reg_volume.match(line)
        self.surface = self._reg_surface.match(line)
        self.sav = self._reg_sav.match(line)
        self.edginess = self._reg_edginess.match(line)
        self.compactess = self._reg_compactess.match(line)
        self.mahalquant = self._reg_mahalquant.match(line)
        self.mahalmin = self._reg_mahalmin.match(line)
        self.mahalmean = self._reg_mahalmean.match(line)
        self.mahalmax = self._reg_mahalmax.match(line)
        self.mean = self._reg_mean.match(line)
        self.variance = self._reg_variance.match(line)
        self.distancegrav = self._reg_distancegrav.match(line)
        self.propneightwm = self._reg_propneightwm.match(line)
        self.propneighgm = self._reg_propneighgm.match(line)
        self.propneighcsf = self._reg_propneighcsf.match(line)
        self.propneighnb = self._reg_propneighnb.match(line)
        self.propneighout = self._reg_propneighout.match(line)
        self.grav = self._reg_grav.match(line)
        self.diffgrav = self._reg_diffgrav.match(line)
        self.dgm = self._reg_dgm.match(line)
        self.sppot = self._reg_sppot.match(line)
        self.propicsf = self._reg_propicsf.match(line)
        self.propecsf = self._reg_propecsf.match(line)
        self.propcgm = self._reg_propcgm.match(line)
        self.propdgm = self._reg_propdgm.match(line)
        self.propwmh = self._reg_propwmh.match(line)
        self.propartefact = self._reg_propartefact.match(line)
        self.propneighartefact = self._reg_propneighartefact.match(line)
        self.propborderextent = self._reg_propborderextent.match(line)
        self.distventr = self._reg_distventr.match(line)
        self.distgmi = self._reg_distgmi.match(line)
        self.distwmi = self._reg_distwmi.match(line)
        self.distcsfi = self._reg_distcsfi.match(line)
        self.distouti = self._reg_distouti.match(line)
        self.distsp = self._reg_distsp.match(line)
        self.ratiobb = self._reg_ratiobb.match(line)
        self.bbrelgrav = self._reg_bbrelgrav.match(line)
        self.extentbox = self._reg_extentbox.match(line)
        self.ratioext = self._reg_ratioext.match(line)
        self.eigenvalues = self._reg_eigenvalues.match(line)
        self.ratioeigen = self._reg_ratioeigen.match(line)
        self.proporigin = self._reg_proporigin.match(line)
        self.texture = self._reg_texture.match(line)
        self.texture1 = self._reg_texture1.match(line)


def main(argv):
    path_save = None
    id = None
    modalities = 2
    parse_file = None

    try:
        opts, args = getopt.getopt(argv, "hf:p:i:m:", ["file=",
                                                       "modalities=", "id=",
                                                       "path_save="])
    except getopt.GetoptError:
        print('test.py -a <association_file> -p <pattern> ')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('test.py -i <regexp> -m <regexp> ')
            sys.exit()
        elif opt in ("-f", "--file"):
            parse_file = arg
        elif opt in ("-m," ,"--modalities"):
            modalities = arg
        elif opt in ("-p", "--path_save"):
            path_save = arg
        elif opt in ("-i", "--id"):
            id = arg

    data = parse(parse_file, modalities)
    if path_save is None:
        path_save = os.path.split(parse_file)[0]
    if id is None:
        id = os.path.split(parse_file)[1]
        id = os.path.splitext(id)[0]
    file_save = os.path.join(path_save, 'ParsedLesion_'+id+'.csv')
    data = pd.DataFrame(data)
    data.to_csv(file_save)


if __name__ == '__main__':
    main(sys.argv[1:])





