import hashlib
import shutil
import uuid as uuidlib
from pathlib import Path
from xml.sax.saxutils import quoteattr

import numpy as np

MZ_DTYPE = np.dtype("<f8")
INTENSITY_DTYPE = np.dtype("<f4")


def _cv(ref: str, acc: str, name: str, value: str = "", unit: tuple[str, str] | None = None) -> str:
    out = f'<cvParam cvRef="{ref}" accession="{acc}" name={quoteattr(name)} value={quoteattr(value)}'
    if unit:
        out += f' unitCvRef="UO" unitAccession="{unit[0]}" unitName="{unit[1]}"'
    return out + "/>"


class ImzMLWriter:
    def __init__(
        self,
        path: str | Path,
        mode: str = "processed",
        mz: np.ndarray | None = None,
        pixel_size: tuple[float, float] | None = None,
        centroided: bool = True,
        source_name: str = "open_slx",
    ):
        if mode not in ("continuous", "processed"):
            raise ValueError("mode must be continuous or processed")
        if mode == "continuous" and mz is None:
            raise ValueError("continuous mode needs a shared m/z axis")
        path = Path(path)
        self.imzml_path = path.with_name(path.name + ".imzML")
        self.ibd_path = path.with_name(path.name + ".ibd")
        self._body_path = path.with_name(path.name + ".imzML.body.tmp")
        self.mode = mode
        self.pixel_size = pixel_size
        self.centroided = centroided
        self.source_name = source_name
        self.uuid = uuidlib.uuid4()
        self.sha1 = hashlib.sha1()
        self._ibd = open(self.ibd_path, "wb")
        try:
            self._body = open(self._body_path, "w", encoding="utf-8")
        except BaseException:
            self._ibd.close()
            self.ibd_path.unlink(missing_ok=True)
            raise
        self._pos = 0
        self._write(self.uuid.bytes)
        self.n_spectra = 0
        self._shared: tuple[int, int] | None = None
        if mode == "continuous":
            self._n_shared = len(mz)
            data = np.asarray(mz, MZ_DTYPE).tobytes()
            self._shared = (self._pos, len(mz))
            self._write(data)
        self._max_x = 0
        self._max_y = 0
        self._closed = False

    def _write(self, data: bytes) -> None:
        self._ibd.write(data)
        self.sha1.update(data)
        self._pos += len(data)

    def add_spectrum(self, x: int, y: int, mz: np.ndarray | None, intensities: np.ndarray) -> None:
        inten = np.asarray(intensities, INTENSITY_DTYPE)
        if self.mode == "continuous":
            if len(inten) != self._n_shared:
                raise ValueError("intensity length differs from the shared m/z axis")
            mz_off, mz_len = self._shared
        else:
            mz = np.asarray(mz, MZ_DTYPE)
            if len(mz) != len(inten):
                raise ValueError("m/z and intensity lengths differ")
            mz_off, mz_len = self._pos, len(mz)
            self._write(mz.tobytes())
        int_off = self._pos
        self._write(inten.tobytes())
        self._body.write(self._spectrum_xml(self.n_spectra, x, y, mz_off, mz_len, int_off, len(inten)))
        self.n_spectra += 1
        self._max_x = max(self._max_x, x)
        self._max_y = max(self._max_y, y)

    def close(self) -> None:
        if self._closed:
            return
        self._closed = True
        self._ibd.close()
        self._body.close()
        with open(self.imzml_path, "w", encoding="utf-8") as out:
            out.write(self._header_xml())
            with open(self._body_path, encoding="utf-8") as body:
                shutil.copyfileobj(body, out)
            out.write("</spectrumList></run></mzML>\n")
        self._body_path.unlink()

    def discard(self) -> None:
        self._closed = True
        self._ibd.close()
        self._body.close()
        for p in (self.ibd_path, self._body_path, self.imzml_path):
            p.unlink(missing_ok=True)

    def __enter__(self) -> "ImzMLWriter":
        return self

    def __exit__(self, *exc) -> None:
        self.close()

    def _spectrum_xml(self, i: int, x: int, y: int, mo: int, ml: int, io: int, il: int) -> str:
        lines = [
            f'<spectrum id="Scan={i + 1}" defaultArrayLength="0" index="{i}">',
            '<cvParam cvRef="MS" accession="MS:1000579" name="MS1 spectrum" value=""/>',
            '<scanList count="1">',
            '<cvParam cvRef="MS" accession="MS:1000795" name="no combination" value=""/>',
            "<scan>",
            _cv("IMS", "IMS:1000050", "position x", str(x)),
            _cv("IMS", "IMS:1000051", "position y", str(y)),
            "</scan></scanList>",
            '<binaryDataArrayList count="2">',
            '<binaryDataArray encodedLength="0"><referenceableParamGroupRef ref="mzArray"/>',
            _cv("IMS", "IMS:1000103", "external array length", str(ml)),
            _cv("IMS", "IMS:1000104", "external encoded length", str(ml * MZ_DTYPE.itemsize)),
            _cv("IMS", "IMS:1000102", "external offset", str(mo)),
            "<binary/></binaryDataArray>",
            '<binaryDataArray encodedLength="0"><referenceableParamGroupRef ref="intensities"/>',
            _cv("IMS", "IMS:1000103", "external array length", str(il)),
            _cv("IMS", "IMS:1000104", "external encoded length", str(il * INTENSITY_DTYPE.itemsize)),
            _cv("IMS", "IMS:1000102", "external offset", str(io)),
            "<binary/></binaryDataArray>",
            "</binaryDataArrayList></spectrum>",
        ]
        return "\n".join(lines) + "\n"

    def _header_xml(self) -> str:
        ps = self.pixel_size
        lines = [
            '<?xml version="1.0" encoding="UTF-8"?>',
            '<mzML xmlns="http://psi.hupo.org/ms/mzml" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" '
            'xsi:schemaLocation="http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.0.xsd" version="1.1">',
            '<cvList count="3">',
            '<cv id="MS" fullName="Proteomics Standards Initiative Mass Spectrometry Ontology" version="4.1.0" URI="https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo"/>',
            '<cv id="UO" fullName="Unit Ontology" version="1.0" URI="https://raw.githubusercontent.com/bio-ontology-research-group/unit-ontology/master/unit.obo"/>',
            '<cv id="IMS" fullName="Imaging MS Ontology" version="1.1.0" URI="https://raw.githubusercontent.com/imzML/imzML/master/imagingMS.obo"/>',
            "</cvList>",
            '<fileDescription><fileContent>',
            _cv("MS", "MS:1000579", "MS1 spectrum"),
            _cv("MS", "MS:1000127" if self.centroided else "MS:1000128", "centroid spectrum" if self.centroided else "profile spectrum"),
            _cv("IMS", "IMS:1000080", "universally unique identifier", "{" + str(self.uuid).upper() + "}"),
            _cv("IMS", "IMS:1000091", "ibd SHA-1", self.sha1.hexdigest().upper()),
            _cv("IMS", "IMS:1000030", "continuous") if self.mode == "continuous" else _cv("IMS", "IMS:1000031", "processed"),
            "</fileContent></fileDescription>",
            '<referenceableParamGroupList count="2">',
            '<referenceableParamGroup id="mzArray">',
            _cv("MS", "MS:1000523", "64-bit float"),
            _cv("MS", "MS:1000576", "no compression"),
            _cv("MS", "MS:1000514", "m/z array"),
            _cv("IMS", "IMS:1000101", "external data", "true"),
            "</referenceableParamGroup>",
            '<referenceableParamGroup id="intensities">',
            _cv("MS", "MS:1000521", "32-bit float"),
            _cv("MS", "MS:1000576", "no compression"),
            _cv("MS", "MS:1000515", "intensity array"),
            _cv("IMS", "IMS:1000101", "external data", "true"),
            "</referenceableParamGroup>",
            "</referenceableParamGroupList>",
            '<softwareList count="1">',
            f'<software id="open_slx" version="0.1.0">{_cv("MS", "MS:1000799", "custom unreleased software tool", self.source_name)}</software>',
            "</softwareList>",
            '<scanSettingsList count="1"><scanSettings id="scansettings1">',
            _cv("IMS", "IMS:1000042", "max count of pixels x", str(self._max_x)),
            _cv("IMS", "IMS:1000043", "max count of pixels y", str(self._max_y)),
        ]
        if ps:
            lines += [
                _cv("IMS", "IMS:1000044", "max dimension x", f"{round(self._max_x * ps[0])}", ("UO:0000017", "micrometer")),
                _cv("IMS", "IMS:1000045", "max dimension y", f"{round(self._max_y * ps[1])}", ("UO:0000017", "micrometer")),
                _cv("IMS", "IMS:1000046", "pixel size (x)", f"{ps[0]:g}", ("UO:0000017", "micrometer")),
                _cv("IMS", "IMS:1000047", "pixel size y", f"{ps[1]:g}", ("UO:0000017", "micrometer")),
            ]
        lines += [
            "</scanSettings></scanSettingsList>",
            '<instrumentConfigurationList count="1"><instrumentConfiguration id="IC1">',
            _cv("MS", "MS:1000031", "instrument model"),
            "</instrumentConfiguration></instrumentConfigurationList>",
            '<dataProcessingList count="1"><dataProcessing id="export">',
            '<processingMethod order="0" softwareRef="open_slx">' + _cv("MS", "MS:1000544", "Conversion to mzML") + "</processingMethod>",
            "</dataProcessing></dataProcessingList>",
            f'<run id="run1" defaultInstrumentConfigurationRef="IC1"><spectrumList count="{self.n_spectra}" defaultDataProcessingRef="export">',
        ]
        return "\n".join(l for l in lines if l) + "\n"
