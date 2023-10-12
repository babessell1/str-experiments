import os
import boto3
import tarfile
import pandas as pd
import glob
import csv


def download_tar_files_from_directory(bucket_name, directory_name, destination_folder, unpack_directory="output"):
    s3 = boto3.client('s3')

    # List all objects in the directory
    response = s3.list_objects_v2(Bucket=bucket_name, Prefix=directory_name)

    # Iterate over the objects and download tar files
    # paginator
    paginator = s3.get_paginator('list_objects_v2')
    for response in paginator.paginate(Bucket=bucket_name, Prefix=directory_name):
        for obj in response.get('Contents', []):
            key = obj['Key']
            filename = os.path.basename(key)
            samples = filename.split('.')[0].split('___')

            # Check if the object is a tar file
            if filename.endswith('.tar'):
                if not all([len(glob.glob(os.path.join(destination_folder, unpack_directory, f'{sample}*'))) > 0 for sample in samples]):
                    # Download the tar file
                    destination_path = os.path.join(destination_folder, filename)
                    print(f'Downloading {filename}', end='\r')
                s3.download_file(bucket_name, key, destination_path)
            

def unpack_tar_files(directory, unpack_directory="output"):
    for filename in os.listdir(directory):
        file_path = os.path.join(directory, filename)

        # only unpack if it's a tar file and not already unpacked
        if filename.endswith('.tar'):
            samples = filename.split('.')[0].split('___')
            # check that not all samples are already unpacked at output/{sample}*
            if not all([len(glob.glob(os.path.join(directory, unpack_directory, f'{sample}*'))) > 0 for sample in samples]):
                print(f'Unpacking {filename}')
                with tarfile.open(file_path, 'r') as tar:
                    try:
                        tar.extractall(directory)
                    except:
                        print(f'Error unpacking {filename}')
                        print("--------------------------------")
                    print(f'Unpacked {filename}', end='\r')
                os.remove(file_path)