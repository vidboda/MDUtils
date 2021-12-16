import os
import urllib3
from MobiDetailsApp import config, md_utilities

# script to test VV API and get current URL for VV
# david 20210217


main():
    http = urllib3.PoolManager(cert_reqs='CERT_REQUIRED', ca_certs=certifi.where())
    try:
        hello = json.loads(http.request('GET', https://rest.variantvalidator.org/hello/?content-type=application/json).data.decode('utf-8')),
        if hello['status'] == "hello_world":
            with open("test.txt", "w") as vv_url_file:
                fo.write("This is Test Data")

if __name__ == '__main__':
    main()
