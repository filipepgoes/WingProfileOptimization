import subprocess

def clmax(commands_file: str): 
    subprocess.run(["process.exe", "starter.exe", "20", "starter.exe",commands_file],cwd="xfoil",shell=True)

def main():
    clmax("entrada.key")

if __name__=='__main__':
    main()