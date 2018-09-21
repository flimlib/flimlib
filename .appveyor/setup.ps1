If ($Env:PLATFORM -match "x64") {
    $Env:MINGW64_BIN="C:\mingw-w64\x86_64-6.3.0-posix-seh-rt_v5-rev1\mingw64\bin"
    $Env:MINGW32_BIN="C:\mingw-w64\x86_64-6.3.0-posix-seh-rt_v5-rev1\x86_64-w64-mingw32\bin"
} ElseIf ($Env:PLATFORM -match "x86") {
    $Env:MINGW64_BIN="C:\mingw-w64\i686-6.3.0-posix-dwarf-rt_v5-rev1\mingw32\bin"
    $Env:MINGW32_BIN="C:\mingw-w64\i686-6.3.0-posix-dwarf-rt_v5-rev1\mingw32\i686-w64-mingw32\bin"
}

# Clean up the machine-wide Path environment variable to a minimum (to enable GPG4Win installation)
[Environment]::SetEnvironmentVariable("Path", "C:\Program Files\AppVeyor\BuildAgent\", "Machine");

# Set up this process' Path environment variable
$Env:Path = "$Env:Path;$Env:GIT_BIN;$Env:MAVEN_BIN;${Env:JAVA_HOME}bin;$Env:GNUPG_HOME;$Env:MINGW64_BIN;$Env:MINGW32_BIN;C:\Windows\system32;C:\Windows;C:\Windows\System32\Wbem;"

# Set up the machine-wide Path environment variable
[Environment]::SetEnvironmentVariable("Path", $Env:Path, "Machine");

"== Path (Execution Environment) =="
$Env:Path

"`n== Path (Machine) =="
[Environment]::GetEnvironmentVariable("Path", "Machine")

"`n== Maven =="
& "mvn" "-B" "--version"

"`n== Java =="
& "cmd.exe" /c 'java -version 2>&1' # Output goes to STDERR which makes this invocation fail the build in PS

# Install GPG4Win
# & ((Split-Path $MyInvocation.InvocationName) + "\install-gpg4win.ps1")

"`n== GPG4Win =="
# & "gpg2.exe" "--version"

# Import the GPG signing key.
$keyFile = "signingkey.asc"
If ((Test-Path Env:\CERTIFICATE_KEY) -and !(Test-Path Env:\APPVEYOR_PULL_REQUEST_NUMBER) -and (Test-Path ("$Env:APPVEYOR_BUILD_FOLDER\.appveyor\$keyFile"))) {
    "`n== Importing GPG keypair =="
    & "cmd.exe" /c "nuget install secure-file -ExcludeVersion 2>&1"
    & "cmd.exe" /c "secure-file\tools\secure-file -decrypt $Env:APPVEYOR_BUILD_FOLDER\.appveyor\signingkey.asc.enc -secret $Env:CERTIFICATE_KEY 2>&1"
    & "cmd.exe" /c "gpg2.exe --batch --fast-import $Env:APPVEYOR_BUILD_FOLDER\.appveyor\$keyFile 2>&1"
}

# Populate the settings.xml configuration.
If (!(Test-Path "$Env:USERPROFILE\.m2")) {
    # Create .m2 directory
    "`n== Creating .m2 directory in user folder =="
    New-Item -Path "$Env:USERPROFILE" -Name ".m2" -ItemType "directory"
}
$settingsFile="$Env:USERPROFILE\.m2\settings.xml"
$customSettings=".appveyor\settings.xml"
If (Test-Path $customSettings) {
    # Copy custom settings.xml
    "`n== Copying .appveyor\settings.xml =="
    Move-Item -Path $customSettings -Destination $settingsFile
} Else {
    # Create default settings.xml
    "`n== Using default settings.xml =="
    $defaultSettings='<settings>
  <servers>
    <server>
      <id>imagej.releases</id>
      <username>appveyor</username>
      <password>${env.MAVEN_PASS}</password>
    </server>
    <server>
      <id>imagej.snapshots</id>
      <username>appveyor</username>
      <password>${env.MAVEN_PASS}</password>
    </server>
  </servers>
  <profiles>
    <profile>
      <id>gpg</id>
      <activation>
        <file>
          <exists>${env.USERPROFILE}\AppData\Roaming\gnupg\pubring.gpg</exists>
        </file>
      </activation>
      <properties>
        <gpg.keyname>${env.GPG_KEY_NAME}</gpg.keyname>
        <gpg.passphrase>${env.GPG_PASSPHRASE}</gpg.passphrase>
      </properties>
    </profile>
  </profiles>
  <activeProfiles>
    <activeProfile>copy-launchers</activeProfile>
  </activeProfiles>
</settings>'
    $defaultSettings | Out-File $settingsFile
}
