import matlab.net.*
import matlab.net.http.*

r = RequestMessage;
uri = URI('https://celestrak.org/NORAD/elements/gp.php?GROUP=iridium&FORMAT=json');
resp = send(r,uri);
status = resp.StatusCode;
data = resp.Body.Data;

data(1)