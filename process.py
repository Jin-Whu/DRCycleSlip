#!/usr/bin/env python3
# coding:utf-8

import sys
import os
import copy
import datetime
import numpy
import matplotlib.pyplot as plt


class Obstype(object):
    """Observation type.

    Attributes:
        interval:epoch interval
    """

    @property
    def epoch(self):
        return self._time

    @epoch.setter
    def epoch(self, value):
        t = [int(float(x)) for x in value.split()]
        self._time = datetime.datetime(t[0], t[1], t[2], t[3], t[4], t[5])

    @property
    def epochnum(self):
        return self._epochnum

    @epochnum.setter
    def epochnum(self, value):
        self._epochnum = value

    interval = None


class GPS(Obstype):
    """GPS observation."""
    C1 = 'C1'
    L1 = 'L1'
    C2 = 'C2'
    L2 = 'L2'
    C5 = 'C5'
    L5 = 'L5'
    f1 = 1575.420 * (10**6)
    f2 = 1227.600 * (10**6)
    f3 = 1176.450 * (10**6)
    w1 = 299792458 / f1
    w2 = 299792458 / f2
    w3 = 299792458 / f3

    def setvalue(self, line):
        """Set observation value.

        Arg:
            line:RINEX observation line.
        """
        obs = line[3:]
        self.c1 = float(obs[self.C1 * 16:self.C1 * 16 + 14])
        self.c2 = float(obs[self.C2 * 16:self.C2 * 16 + 14])
        self.c3 = float(obs[self.C5 * 16:self.C5 * 16 + 14])
        self.l1 = float(obs[self.L1 * 16:self.L1 * 16 + 14])
        self.l2 = float(obs[self.L2 * 16:self.L2 * 16 + 14])
        self.l3 = float(obs[self.L5 * 16:self.L5 * 16 + 14])


class BDS(Obstype):
    """BDS observation."""
    C1 = 'C1'
    L1 = 'L1'
    C2 = 'C2'
    L2 = 'L2'
    C6 = 'C6'
    L6 = 'L6'
    C7 = 'C7'
    L7 = 'L7'
    f1 = 1561.098 * (10**6)
    f2 = 1207.140 * (10**6)
    f3 = 1268.520 * (10**6)
    w1 = 299792458 / f1
    w2 = 299792458 / f2
    w3 = 299792458 / f3

    def setvalue(self, line):
        """Set observation value.

        Arg:
            line:RINEX observation line.
        """
        obs = line[3:]
        try:
            self.c1 = float(obs[self.C2 * 16:self.C2 * 16 + 14])
        except:
            self.c1 = float(obs[self.C1 * 16:self.C1 * 16 + 14])
        self.c2 = float(obs[self.C7 * 16:self.C7 * 16 + 14])
        self.c3 = float(obs[self.C6 * 16:self.C6 * 16 + 14])
        try:
            self.l1 = float(obs[self.L2 * 16:self.L2 * 16 + 14])
        except:
            self.l1 = float(obs[self.L1 * 16:self.L1 * 16 + 14])
        self.l2 = float(obs[self.L7 * 16:self.L7 * 16 + 14])
        self.l3 = float(obs[self.L6 * 16:self.L6 * 16 + 14])


class ComObs(object):
    """Combination observation."""

    def __init__(self, obs, state):
        """Initialize combination observation.

        Args:
            obs:observation.
            state:cycle slip detect state.
        """
        if isinstance(obs, GPS):
            if state == 'static':
                self.combinations = [[-6, 1, 7], [3, 0, -4], [4, -8, 3]]
                self.wavelengths = [29.305, 14.653, 29.305]
                self.threshold = [0.744, 0.412, 0.756]
            elif state == 'move':
                self.combinations = [[-6, 1, 7], [3, 0, -4], [4, -8, 3]]
                self.wavelengths = [29.305, 14.653, 29.305]
                self.threshold = [0.744, 0.412, 0.756]
                # self.combinations = [[-6, 1, 7], [3, 0, -4], [-1, 8, -7]]
                # self.wavelengths = [29.305, 14.653, 29.305]
                # self.threshold = [0.744, 0.412, 0.976]
        elif isinstance(obs, BDS):
            if state == 'static':
                self.combinations = [[-4, 1, 4], [-3, 6, -2], [4, -2, -3]]
                self.wavelengths = [8.14, 13.321, 12.211]
                self.threshold = [0.492, 0.568, 0.444]
            elif state == 'move':
                self.combinations = [[-1, -5, 6], [5, 3, -9], [7, -8, -1]]
                self.wavelengths = [20.932, 29.305, 146.526]
                self.threshold = [0.912, 0.98, 0.86]
                # self.combinations = [[-8, 3, 7], [5, 3, -9], [7, -8, -1]]
                # self.wavelengths = [24.421, 29.305, 146.526]
                # self.threshold = [1.052, 0.98, 0.86]

        self.comvalues = list()
        compr = sum([i * j
                     for i, j in zip([1 / 3.0] * 3, [obs.c1, obs.c2, obs.c3])])

        for combination, wavelength in zip(self.combinations,
                                           self.wavelengths):
            comphase = sum(
                [i * j for i, j in zip(combination, [obs.l1, obs.l2, obs.l3])])
            self.comvalues.append(comphase - compr / wavelength)


class Cycleslip(object):
    """Cycleslip slip"""

    def __init__(self, epoch, epochnum, cycleslip, deltaN):
        """Initialize.

        Args:
            epoch:cycle slip epoch.
            epochnum:cycle slip epoch num.
            cycleslip:cycle slip.
            deltaN:deltaN.
        """
        self.epoch = epoch
        self.epochnum = epochnum
        self.cycleslip = cycleslip
        self.deltaN = deltaN

    def __str__(self):
        """Prin cycleslip."""
        return '%s %s %s %s %.2f %.2f %.2f' % (
            self.epochnum, self.cycleslip[0], self.cycleslip[1],
            self.cycleslip[2], *self.deltaN)


class DeltaN(object):
    """Cycle slip detection value."""

    threshold = None
    combinations = None

    def __init__(self, epochnum, value):
        """Initialize cycle slip detection.

        Args:
            epochnum:epoch number.
            value:detection value.
        """
        self.epochnum = epochnum
        self.value = value


class IonDiff(object):
    """Ionsphere errors difference, including first-order and second-order."""

    def __init__(self, epochnum, diff, diff2):
        """Initialize IonDiff.

        Args:
            epochnum:epoch num.
            diff:first-order difference.
            diff2:second-order difference.
        """
        self.epochnum = epochnum
        self.diff = diff
        self.diff2 = diff2


def readheader(f, system):
    """Read observation file header.

    stract observation position from header.

    Arg:
        f:file handle.
        system:GNSS system.
    """
    for line in f:
        if 'SYS / # / OBS TYPES' in line and line.startswith(system):
            obs = line[7:line.index('SYS')].split()
            for i in range(len(obs)):
                if system == 'G':
                    obstype = GPS
                if system == 'C':
                    obstype = BDS

                if getattr(obstype, obs[i][:2], None):
                    setattr(obstype, obs[i][:2], i)

        if 'INTERVAL' in line:
            Obstype.interval = float(line.split()[0])

        if 'END OF HEADER' in line:
            break


def readcycleslip(prn):
    """Read cycle slip from file.

    Arg:
        prn:satellite prn.
    """
    global precycleslip
    with open(
            os.path.join(
                os.path.dirname(__file__), ''.join([prn, 'cycleslip.txt'
                                                    ]))) as f:
        next(f)
        for line in f:
            if line.strip() == '':
                break
            temp = list(map(int, line.split()))
            precycleslip.append(Cycleslip(None, temp[0], temp[1:], None))


def insertcycleslip(obs):
    """Insert cycle slip in observation.

    Arg:
        obs:observation.
    """
    global precycleslip
    for slip in precycleslip:
        if obs.epochnum >= slip.epochnum:
            obs.l1 += slip.cycleslip[0]
            obs.l2 += slip.cycleslip[1]
            obs.l3 += slip.cycleslip[2]


def iondiff(obs3epoch):
    """Calculate ionsphere difference.

    Arg:
        obs3epoch:three successive observations of epoch.
    """
    global iond

    # calculate ionsphere errors difference
    diff = obs3epoch[2].f3**2 * (
        obs3epoch[2].w1 *
        (obs3epoch[2].l1 - obs3epoch[1].l1) - obs3epoch[2].w3 *
        (obs3epoch[2].l3 - obs3epoch[1].l3)) / (
            obs3epoch[2].f1**2 - obs3epoch[2].f3**2)
    diff2 = diff - obs3epoch[2].f3**2 * (
        obs3epoch[2].w1 *
        (obs3epoch[1].l1 - obs3epoch[0].l1) - obs3epoch[2].w3 *
        (obs3epoch[1].l3 - obs3epoch[0].l3)) / (obs3epoch[2].f1**2 -
                                                obs3epoch[2].f3**2)
    iond.append(IonDiff(obs3epoch[2].epochnum, diff, diff2))


def detect(obs3epoch, system, state):
    """Detect cycleslip.

    Args:
        obs3epoch:three successive observations of epoch.
        system:GNSS system.
        state:cycle slip detect state.
    """
    global cycleslip
    global deltaN_lst

    # detect cycle slip
    comobslst = [ComObs(obs, state) for obs in obs3epoch]
    DeltaN.threshold = comobslst[0].threshold
    DeltaN.combinations = comobslst[0].combinations
    deltaN = list()
    hascycleslip = False
    for i in range(3):
        threshold = comobslst[0].threshold[i]
        value = comobslst[2].comvalues[i] - 2 * comobslst[1].comvalues[
            i] + comobslst[0].comvalues[i]
        deltaN.append(value)
        deltaN_lst.append(DeltaN(obs3epoch[2].epochnum, value))
        if abs(deltaN[i]) > threshold:
            hascycleslip = True
    if not hascycleslip:
        return
    A = numpy.array(comobslst[0].combinations)
    L = numpy.round(numpy.array(deltaN))
    x = numpy.linalg.solve(A, L)
    # check if cycle slip
    obs3epoch[2].l1 -= x[0]
    obs3epoch[2].l2 -= x[1]
    obs3epoch[2].l3 -= x[2]
    comobslst = [ComObs(obs, state) for obs in obs3epoch]
    for i in range(3):
        threshold = comobslst[0].threshold[i]
        value = comobslst[2].comvalues[i] - 2 * comobslst[1].comvalues[
            i] + comobslst[0].comvalues[i]
        if abs(value) > threshold:
            hascycleslip = False
    if not hascycleslip:
        obs3epoch[2].l1 += x[0]
        obs3epoch[2].l2 += x[1]
        obs3epoch[2].l3 += x[2]
        return
    cycleslip.append(
        Cycleslip(obs3epoch[2].epoch, obs3epoch[2].epochnum, x, deltaN))


def repair(obs):
    """Repair observation."""
    global cycleslip
    for slip in cycleslip:
        if obs.epoch >= slip.epoch:
            obs.l1 -= slip.cycleslip[0]
            obs.l2 -= slip.cycleslip[1]
            obs.l3 -= slip.cycleslip[2]


def plotcycleslip(prn):
    """Plot cycleslip.

    Arg:
        prn:satellite prn.
    """
    global cycleslip
    global deltaN_lst
    plt.style.use('ggplot')
    f, axes = plt.subplots(3, 1, sharex=True)
    for i in range(3):
        ax = axes[i]
        combination = DeltaN.combinations[i]
        label = ','.join([str(num) for num in combination])
        label = r'$\Delta N_{(%s)}$' % label
        deltaN_c = [deltaN_lst[i + 3 * j]
                    for j in range(0, int(len(deltaN_lst) / 3))]
        ax.plot(
            [temp.epochnum for temp in deltaN_c],
            [temp.value for temp in deltaN_c],
            label=label)
        line1, line2 = ax.plot([deltaN_c[0].epochnum, deltaN_c[-1].epochnum],
                               [DeltaN.threshold[i], DeltaN.threshold[i]],
                               [deltaN_c[0].epochnum, deltaN_c[-1].epochnum],
                               [-DeltaN.threshold[i], -DeltaN.threshold[i]])
        plt.setp(line2, color=plt.getp(line1, 'color'))
        ax.set_xlim(deltaN_c[0].epochnum, deltaN_c[-1].epochnum)
        ax.set_ylim([-2.0, 2.0])
        ax.set_yticks([-2, -DeltaN.threshold[i], DeltaN.threshold[i], 2])
        plt.setp(ax.xaxis.get_ticklabels(), size=10, weight='bold')
        plt.setp(ax.yaxis.get_ticklabels(), size=10, weight='bold')
        ax.legend()
    axes[2].set_xlabel('Epoch/S', size=12, weight='bold')
    axes[1].set_ylabel('CycleSlip detection value', size=12, weight='bold')
    axes[0].set_title(prn, size=15, weight='bold')
    plt.savefig(
        r'C:\Users\jin\OneDrive\graduationproject\program\picture\%s.eps' %
        prn,
        bbox_inches='tight',
        dpi=600)
    # plt.show()
    plt.clf()
    plt.close()


def plotiond(prn):
    """Plot ionsphere difference."""
    global iond
    plt.style.use('ggplot')
    epoch = [d.epochnum for d in iond]
    diff = [d.diff for d in iond]
    diff2 = [d.diff2 for d in iond]
    plt.plot(epoch, diff, label=r'$\Delta I$')
    plt.plot(epoch, diff2, label=r'$\Delta\Delta I$')
    plt.xlim([epoch[0], epoch[-1]])
    plt.ylim([0.2, -0.2])
    plt.xlabel('Epoch/S', size=12, weight='bold')
    plt.ylabel('Ionspheric delay variation[m]', size=12, weight='bold')
    ax = plt.gca()
    plt.setp(ax.xaxis.get_ticklabels(), size=10, weight='bold')
    plt.setp(ax.yaxis.get_ticklabels(), size=10, weight='bold')
    plt.title(prn, size=15, weight='bold')
    plt.legend()
    plt.savefig(
        r'C:\Users\jin\OneDrive\graduationproject\program\picture\%sion.eps' %
        prn,
        bbox_inches='tight',
        dpi=600)
    # plt.show()
    plt.clf()
    plt.close()


def process(filepath, prn, start, end, state):
    """read observation file, detect and repair cycle slip.

    Arg:
        filepath:observation file path.
        prn:satellite prn.
        start:start epoch.
        end:end epoch.
        state:cycle slip detect state.
    """
    try:
        f = open(filepath)
    except FileNotFoundError:
        print('Not find %s' % filepath)
        return
    system = prn[0]
    readheader(f, system)

    # read precycleslip
    readcycleslip(prn)

    # detect cycleslip
    count = 0
    obs3epoch = list()
    obs3ion = list()
    for line in f:
        if '>' in line:
            if system == 'G':
                obs = GPS()
            elif system == 'C':
                obs = BDS()
            obs.epoch = line[2:28]
            count += 1
            obs.epochnum = count

            if count > end:
                break

        if count > start and (
                line.startswith(prn) or
                line.startswith('{0}{1:>2d}'.format(prn[0], int(prn[1:])))):

            try:
                obs.setvalue(line)
            except:
                continue

            obs3ion.append(copy.deepcopy(obs))
            # insert precycleslip
            insertcycleslip(obs)
            obs3epoch.append(obs)

        if len(obs3epoch) == 3:
            # calculate ionsphere difference
            iondiff(obs3ion)
            # detect and repair cycelslip
            repair(obs)
            detect(obs3epoch, system, state)
            obs3epoch.pop(0)
            obs3ion.pop(0)

    #print cycleslip
    global cycleslip
    for slip in cycleslip:
        print(slip)

    # plot
    plotcycleslip(prn)
    plotiond(prn)


if __name__ == '__main__':
    # if len(sys.argv) != 6:
    #     print('\nArgs:')
    #     print('\tfilepath:observation file path.')
    #     print('\tprn:satellite prn.')
    #     print('\tstart:start epoch.')
    #     print('\tend:end epoch.')
    #     print('\tmethod:static or move.')
    #     sys.exit()
    # filepath = sys.argv[1]s
    # prn = sys.argv[2]
    # start = int(sys.argv[3])
    # end = int(sys.argv[4])
    # state = sys.argv[5]
    # process(filepath, prn, start, end, state)
    filepath = [r'C:\Users\jin\Downloads\jfng0760.13o'] * 4
    prns = ['C03', 'C09', 'C12', 'G24']
    startepoch = [1, 1, 1300, 1]
    endepoch = [3000, 1500, 2300, 800]
    states = ['static'] * 4
    for path, prn, start, end, state in zip(filepath, prns, startepoch,
                                            endepoch, states):
        print(prn)
        iond = list()
        cycleslip = list()
        deltaN_lst = list()
        precycleslip = list()
        process(path, prn, start, end, state)
        print()
