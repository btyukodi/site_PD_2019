from flask import Flask, render_template
import data_process as dp
from bson.objectid import ObjectId

app = Flask(__name__)

@app.route('/<run_id>')
def index(run_id):
  dp.ovito(ObjectId(run_id))
  #print "aaa"
  return '<!DOCTYPE html><html><body onload="self.close()"></body></html>'#render_template('template.html')

##@app.route('/my-link/')
##def my_link():
##  print 'I got clicked!'
##
##  return 'Click.'

if __name__ == '__main__':
  app.run(debug=True)
