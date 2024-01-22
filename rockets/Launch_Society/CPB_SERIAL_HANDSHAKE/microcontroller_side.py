import time
import supervisor

while True:
	if supervisor.runtime.serial_bytes_available:
		df = input.strip()
		print('Received df command = ',df)
		#Tell the servo to move
	time.sleep(0.01)