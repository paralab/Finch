import os
import shutil


def FileGenerator(filename, template):
    def generator(case, job_dir):
        path = os.path.join(job_dir, filename)
        with open(path, 'w') as f:
            f.write(template.format(**case))

    return generator


def FileCopyGenerator(destpath, srcpath):
    assert os.path.exists(srcpath), "File copy generator: source file '{}' does not exist".format(srcpath)

    def generator(case, job_dir):
        shutil.copyfile(srcpath, os.path.join(job_dir, destpath))

    return generator


def FileLinkGenerator(destpath, srcpath):
    assert os.path.exists(srcpath), "File link generator: source file '{}' does not exist".format(srcpath)

    def generator(case, job_dir):
        os.symlink(srcpath, os.path.join(job_dir, destpath))

    return generator


def FolderCopyGenerator(destpath, srcpath):
    assert os.path.exists(srcpath), "Folder copy generator: source folder '{}' does not exist".format(srcpath)

    def generator(case, job_dir):
        shutil.copytree(srcpath, os.path.join(job_dir, destpath))

    return generator
