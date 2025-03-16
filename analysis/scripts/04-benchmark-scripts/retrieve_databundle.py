from zipfile import ZipFile
import io
import requests

if __name__ == "__main__":

    url = "https://zenodo.org/records/15036329/files/2023-osier-simulation-results.zip?download=1"
    download_path = "../data/osier_data"

    print("Downloading databundle from Zenodo")

    r = requests.get(url)

    if r.status_code == 200:
        print('Download successful! Unzipping file bundle.')
        z = ZipFile(io.BytesIO(r.content))

        z.extractall("../data/")