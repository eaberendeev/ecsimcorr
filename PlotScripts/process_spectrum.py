import numpy as np
import argparse
import os

def main():
    parser = argparse.ArgumentParser(description='Aggregate particle counts within a specified energy range across multiple files.')
    parser.add_argument('--start', type=int, required=True, help='Starting file number')
    parser.add_argument('--end', type=int, required=True, help='Ending file number (inclusive)')
    parser.add_argument('--emin', type=float, required=True, help='Minimum energy value')
    parser.add_argument('--emax', type=float, required=True, help='Maximum energy value')
    parser.add_argument('--input', type=str, default='energy_spectrum_{:06d}.txt',
                        help='Input filename template (e.g., energy_spectrum_{:06d}.txt)')
    parser.add_argument('--output', type=str, default='output.txt',
                        help='Output filename for time vs count data')
    
    args = parser.parse_args()

    with open(args.output, 'w') as outfile:
        outfile.write('Time ParticleCount\n')
        
        for time in range(args.start, args.end + 1):
            filename = args.input.format(time)
            
            if not os.path.exists(filename):
                print(f"Warning: File {filename} not found. Skipping time point {time}.")
                continue
            
            try:
                data = np.loadtxt(filename)
                energies = data[:, 0]
                counts = data[:, 1]
            except Exception as e:
                print(f"Error reading {filename}: {e}. Skipping.")
                continue
            
            # Используем бинарный поиск для границ
            left = np.searchsorted(energies, args.emin, side='left')
            right = np.searchsorted(energies, args.emax, side='right')
            
            total = np.sum(counts[left:right]) if left < right else 0
            
            outfile.write(f'{time} {total}\n')

if __name__ == '__main__':
    main()