#!/usr/bin/env python3

#import
import argparse

#define functions/classes

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbose', type = bool, help = "Set to True to print status updates. Default True", default = True)
    #add args
    
    args = parser.parse_args()
    return args

def main():
    args = argparser()
    #insert code
    pass

if __name__ == "__main__":
    main()