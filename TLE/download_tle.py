
import requests
import json


url = "https://celestrak.org/NORAD/elements/gp.php?GROUP=iridium&FORMAT=json"
resp = requests.get(url)
data = resp.json()
print(data)
