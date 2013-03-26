import sys

def lib_main(): pass

class Progress():

	def __init__(self,total, verb='full', logfile=sys.stderr):
		import time
		self.total = total
		self.current = 0
		self.type = verb
		self.logfile = logfile
		self.ltime = time.time()
		self.lcurrent = self.current
		self.lpercentage = 0

	def __enter__(self):
		if self.type == 'minimal': self.logfile.write('0%                 50%                 100%\n')

	def update(self):
		import time
		self.current += 1
		self.percentage = int(round(100*float(self.current)/self.total))
		if self.percentage % 5 == 0 and self.percentage != self.lpercentage:
			self.stf=int(round((self.total-self.current)/((self.current-self.lcurrent)/(time.time()-self.ltime))))
			if self.type == 'full': self.logfile.write(
				'#Progress => '+str(self.percentage)+'%, '+
				str( round((self.current-self.lcurrent)/(time.time()-self.ltime),2) )+' reads/second, '+
				time.strftime("%A, %d %b %Y %H:%M:%S",time.localtime())+
				', left:'+str(self.stf/60)+'min '+str(self.stf%60)+'s'+
				'\n'
				)
			if self.type == 'minimal': self.logfile.write('..')
			self.ltime = time.time()
			self.lcurrent = self.current
			self.lpercentage = self.percentage

	def __exit__(self, *args):
		self.logfile.write('\n')

if __name__ == "__main__":
    lib_main()
