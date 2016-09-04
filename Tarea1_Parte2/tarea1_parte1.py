from astropy.coordinates import SkyCoord
from astropy.io import ascii


class Tarea1Parte1:

    def __init__(self):
        # Creamos los atributos.
        # Zero points
        self.f555_zp = 25.255
        self.f814_zp = 24.849

        # Reddening
        self.f555_av = 0.103
        self.f814_av = 0.056

        # Magnitudes
        self.f555_mag, self.f814_mag = self._extract()

    def _extract(self):
        names = ['cat_f555.cat', 'cat_f814.cat']  # Nombres de los catálogos
        skies = []  # Lista de SkyCoord.
        cats = []  # Lista de catálogos.

        for name in names:
            data = ascii.read(name)
            # Filtramos los valores 99 de magnitud para objetos muy débiles.
            mask = data['MAG_APER'] != 99.0
            data = data[mask]

            cats.append(data)
            skies.append(SkyCoord(ra=data['ALPHA_J2000'], dec=data[
                         'DELTA_J2000'], unit='deg'))

        # Obtenemos los matches. Se busca el más cercano de cada estrella en
        # skies[0] en skies[1].
        idx, d2d, d3d = SkyCoord.match_to_catalog_sky(skies[0], skies[1])
        # Creamos una lista de tuplas que relaciona ambos catálogos.
        matches = list(zip(range(0, len(idx)), idx, d2d))
        # Quitamos los que están repetidos y tiene una distancia mayor a 0.1
        # arcosegundos:
        temp_matches = []
        seconds = set()
        matches.sort(key=lambda x: x[2], reverse=True)
        for x in matches:
            if x[1] not in seconds and x[2].arcsec < 0.1:
                seconds.add(x[1])
                temp_matches.append(x)
        matches = temp_matches

        # Ahora extraemos de los catálogos.
        idx1, idx2, d2d = zip(*matches)
        cats[0] = cats[0][list(idx1)]
        cats[1] = cats[1][list(idx2)]

        # f555 = list(map(lambda x: -2.5*np.log10(x), cats[0]['FLUX_APER']))
        # f814 = list(map(lambda x: -2.5*np.log10(x), cats[1]['FLUX_APER']))
        f555 = list(map(lambda x: x, cats[0]['MAG_APER']))
        f814 = list(map(lambda x: x, cats[1]['MAG_APER']))

        # Corrección para F555w y F814w
        mag = [f555, f814]
        # m_0 = [25.255, 24.849]
        m_0 = [0, 0]
        apcor = [0.241, 0.425]
        a_v = [0.103, 0.056]

        for idx in [0, 1]:
            mag[idx] = list(map(lambda x: x + m_0[idx] -
                                apcor[idx] - a_v[idx], mag[idx]))
        return (mag[0], mag[1])
