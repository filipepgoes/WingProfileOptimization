import subproces

def set_line(filename: str, line_description: str, line_values: str):
    line_number=None
    with open(filename) as file:
        data = file.readlines()
        for i in range(len(data)):
            line = data[i].rstrip()  
            if line_description == line:
                line_number=i
    if not line_number:
        raise('Description not found.')
    data[line_number+1] = line_value
    with open(filename, 'w') as file:
        file.writelines(data)

def run(): 
    subprocess.run(["blwf28.exe", "aircraft.dat"],cwd="blwf",shell=True)

def main():
    run()

if __name__=='__main__':
    main()