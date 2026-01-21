# C# Naming Conventions Examples

> **Domain: Gaming Platform**
> Examples use Player, Game, Level, Score, Achievement, and Match entities.

## Table of Contents

1. [Variables](#variables)
2. [Methods](#methods)
3. [Classes](#classes)
4. [Interfaces](#interfaces)
5. [Constants and Enums](#constants-and-enums)
6. [Booleans](#booleans)
7. [Collections](#collections)
8. [Magic Numbers](#magic-numbers)

---

## Variables

### Use Intention-Revealing Names

**❌ BAD:**
```csharp
int d; // elapsed time in seconds
int s; // score
string n;
var list = new List<int>();
```

**✅ GOOD:**
```csharp
int elapsedTimeInSeconds;
int playerScore;
string playerName;
var activeGameSessions = new List<GameSession>();
```

### Avoid Single Letter Names

**❌ BAD:**
```csharp
public void Calculate(int a, int b, int c)
{
    var x = a + b;
    var y = x * c;
    return y;
}
```

**✅ GOOD:**
```csharp
public int CalculateFinalScore(int baseScore, int bonusMultiplier, int comboCount)
{
    var adjustedScore = baseScore + bonusMultiplier;
    var finalScore = adjustedScore * comboCount;
    return finalScore;
}
```

**Exception - Loop Variables:**
```csharp
// OK for simple loops
for (int i = 0; i < count; i++) { }

// Better for complex loops
for (int levelIndex = 0; levelIndex < levels.Count; levelIndex++)
{
    var level = levels[levelIndex];
    // ...
}
```

### Avoid Abbreviations

**❌ BAD:**
```csharp
string plyrNm;
int lvlCnt;
decimal xpAmt;
var mgr = new GameManager();
```

**✅ GOOD:**
```csharp
string playerName;
int levelCount;
decimal experiencePoints;
var gameManager = new GameManager();
```

**Exceptions - Well-Known Abbreviations:**
```csharp
// OK - widely understood in gaming
int hp = 100;        // Health Points
int xp = 500;        // Experience Points
int fps = 60;        // Frames Per Second
var npc = new NonPlayerCharacter();
```

### Avoid Hungarian Notation

**❌ BAD:**
```csharp
string strPlayerName;
int iScore;
bool bIsAlive;
decimal decDamage;
List<Player> lstPlayers;
```

**✅ GOOD:**
```csharp
string playerName;
int score;
bool isAlive;
decimal damage;
List<Player> players;
```

### Use Searchable Names

**❌ BAD:**
```csharp
if (status == 2) { }  // What is 2?

var result = value * 60;  // Why 60?
```

**✅ GOOD:**
```csharp
const int PlayerStatusPlaying = 2;
if (status == PlayerStatusPlaying) { }

const int SecondsPerMinute = 60;
var totalSeconds = minutes * SecondsPerMinute;
```

### Avoid Mental Mapping

**❌ BAD:**
```csharp
foreach (var p in players)
{
    var s = p.Score;
    var l = s / 1000;
    // ...
}
```

**✅ GOOD:**
```csharp
foreach (var player in players)
{
    var score = player.Score;
    var level = score / PointsPerLevel;
    // ...
}
```

---

## Methods

### Use Verb Phrases

**❌ BAD:**
```csharp
public void Player() { }
public string Name() { }
public bool Active() { }
```

**✅ GOOD:**
```csharp
public void SpawnPlayer() { }
public string GetPlayerName() { }
public bool IsPlayerActive() { }
```

### Be Specific

**❌ BAD:**
```csharp
public void Process() { }
public void Handle() { }
public void Do() { }
public void Manage() { }
```

**✅ GOOD:**
```csharp
public void ProcessMatchResult() { }
public void HandlePlayerDeath() { }
public void ApplyDamageToEnemy() { }
public void CalculateLevelProgress() { }
```

### Use Consistent Vocabulary

**❌ BAD:**
```csharp
GetPlayer();
FetchGameSession();
RetrieveLeaderboard();
LoadAchievements();
```

**✅ GOOD:**
```csharp
GetPlayer();
GetGameSession();
GetLeaderboard();
GetAchievements();
```

### Command vs Query Names

**Commands (modify state):**
```csharp
public void SavePlayerProgress(Player player) { }
public void UpdateMatchScore(MatchId id, int score) { }
public void DeleteSaveGame(SaveGameId id) { }
public void AwardAchievement(PlayerId playerId, Achievement achievement) { }
```

**Queries (return data):**
```csharp
public Player GetPlayer(PlayerId id) { }
public List<Match> GetActiveMatches() { }
public bool PlayerHasAchievement(PlayerId id, AchievementType type) { }
public int CalculatePlayerRank(PlayerId id) { }
```

### Async Method Names

**✅ GOOD:**
```csharp
public async Task<Player> GetPlayerAsync(PlayerId id) { }
public async Task SaveGameStateAsync(GameState state) { }
public async Task<List<LeaderboardEntry>> GetLeaderboardAsync(int topN) { }
```

---

## Classes

### Use Nouns or Noun Phrases

**❌ BAD:**
```csharp
public class Manage { }
public class Process { }
public class Handle { }
```

**✅ GOOD:**
```csharp
public class GameManager { }
public class ScoreCalculator { }
public class MatchHandler { }
public class AchievementTracker { }
```

### Avoid Generic Names

**❌ BAD:**
```csharp
public class Manager { }
public class Handler { }
public class Helper { }
public class Data { }
```

**✅ GOOD:**
```csharp
public class PlayerRepository { }
public class DamageCalculator { }
public class LeaderboardFormatter { }
public class MatchStatistics { }
```

### Use Design Pattern Names

**✅ GOOD:**
```csharp
public class PlayerFactory { }
public class CharacterBuilder { }
public class AttackStrategy { }
public class PowerUpDecorator { }
public class MatchObserver { }
public class InputAdapter { }
```

---

## Interfaces

### Use 'I' Prefix (C# Convention)

**✅ GOOD:**
```csharp
public interface IPlayerRepository { }
public interface IMatchmakingService { }
public interface ILeaderboardProvider { }
public interface IGameLogger { }
```

### Describe Capability

**✅ GOOD:**
```csharp
public interface IDamageable { }
public interface IMoveable { }
public interface IAttackable { }
public interface ICollectable { }
public interface ISpawnable { }
```

---

## Constants and Enums

### Constants in PascalCase

**❌ BAD:**
```csharp
const int MAX_PLAYERS = 4;
const string DEFAULT_GAME_MODE = "deathmatch";
```

**✅ GOOD:**
```csharp
public const int MaxPlayersPerMatch = 4;
public const string DefaultGameMode = "deathmatch";
```

### Enums and Values

**✅ GOOD:**
```csharp
public enum PlayerStatus
{
    Idle,
    Playing,
    InLobby,
    Spectating,
    Disconnected
}

public enum GameMode
{
    SinglePlayer,
    Cooperative,
    Competitive,
    BattleRoyale
}

public enum DamageType
{
    Physical,
    Fire,
    Ice,
    Lightning,
    Poison
}
```

---

## Booleans

### Use Predicate Phrases

**❌ BAD:**
```csharp
bool alive;
bool flag;
bool status;
```

**✅ GOOD:**
```csharp
bool isAlive;
bool hasShield;
bool canAttack;
bool shouldRespawn;
bool wasDefeated;
```

### Question Form for Methods

**✅ GOOD:**
```csharp
bool IsPlayerEligibleForMatch();
bool HasUnlockedAchievement(AchievementType type);
bool CanPlayerAffordItem(Player player, Item item);
bool ShouldShowTutorial();
bool WasMatchCompleted();
```

### Avoid Negatives

**❌ BAD:**
```csharp
bool isNotAlive;
bool isNotInvincible;

if (!isNotAlive) // Double negative!
```

**✅ GOOD:**
```csharp
bool isAlive;
bool isVulnerable;

if (isAlive) { }
```

---

## Collections

### Use Plural Names

**❌ BAD:**
```csharp
var player = new List<Player>();
var achievement = new List<Achievement>();
```

**✅ GOOD:**
```csharp
var players = new List<Player>();
var achievements = new List<Achievement>();
var levelsById = new Dictionary<int, Level>();
```

### Be Specific About Collection Type

**✅ GOOD:**
```csharp
var activePlayers = new List<Player>();
var pendingMatchRequests = new Queue<MatchRequest>();
var playersByGuildId = new Dictionary<GuildId, List<Player>>();
var unlockedAchievementIds = new HashSet<AchievementId>();
```

---

## Magic Numbers

### Extract to Named Constants

**❌ BAD:**
```csharp
public class MatchValidator
{
    public bool ValidateMatch(Match match)
    {
        if (match.Players.Count > 100)
            return false;

        if (match.Duration > 3600)
            return false;

        if (match.MinLevel < 5)
            return false;

        return true;
    }
}
```

**✅ GOOD:**
```csharp
public class MatchValidator
{
    private const int MaxPlayersPerMatch = 100;
    private const int MaxMatchDurationSeconds = 3600;
    private const int MinimumPlayerLevel = 5;

    public bool ValidateMatch(Match match)
    {
        if (match.Players.Count > MaxPlayersPerMatch)
            return false;

        if (match.Duration > MaxMatchDurationSeconds)
            return false;

        if (match.MinLevel < MinimumPlayerLevel)
            return false;

        return true;
    }
}
```

### Use Configuration

**✅ EVEN BETTER:**
```csharp
public class MatchValidationRules
{
    public int MaxPlayersPerMatch { get; set; } = 100;
    public int MaxMatchDurationSeconds { get; set; } = 3600;
    public int MinimumPlayerLevel { get; set; } = 5;
}

public class MatchValidator
{
    private readonly MatchValidationRules _rules;

    public MatchValidator(MatchValidationRules rules)
    {
        _rules = rules;
    }

    public bool ValidateMatch(Match match)
    {
        return match.Players.Count <= _rules.MaxPlayersPerMatch
            && match.Duration <= _rules.MaxMatchDurationSeconds
            && match.MinLevel >= _rules.MinimumPlayerLevel;
    }
}
```

---

## Real-World Example: Before and After

### Before (Bad Naming)

```csharp
public class GM
{
    private List<P> pl;
    private int s;

    public void U(P p)
    {
        if (p.st == 1)
        {
            var d = 0;
            foreach (var e in p.eq)
            {
                d += e.dmg * e.lvl;
            }

            if (p.vip)
                d = (int)(d * 1.5);

            p.td = d;
            pl.Add(p);
            s += d;
        }
    }
}
```

### After (Good Naming)

```csharp
public class GameManager
{
    private readonly List<Player> _activePlayers;
    private int _totalDamageDealt;

    public void UpdatePlayerStats(Player player)
    {
        if (player.Status == PlayerStatus.Playing)
        {
            var totalDamage = CalculateTotalDamage(player);

            if (player.IsPremiumMember)
                totalDamage = ApplyPremiumBonus(totalDamage);

            player.TotalDamage = totalDamage;
            _activePlayers.Add(player);
            _totalDamageDealt += totalDamage;
        }
    }

    private int CalculateTotalDamage(Player player)
    {
        int totalDamage = 0;
        foreach (var equipment in player.Equipment)
        {
            totalDamage += equipment.Damage * equipment.Level;
        }
        return totalDamage;
    }

    private int ApplyPremiumBonus(int damage)
    {
        const decimal PremiumDamageMultiplier = 1.5m;
        return (int)(damage * PremiumDamageMultiplier);
    }
}
```

---

## Summary

| Principle | Description |
|-----------|-------------|
| **Reveal Intent** | Name explains why it exists and what it does |
| **Avoid Disinformation** | Don't use names that mean something else |
| **Make Distinctions** | Different names mean different things |
| **Use Pronounceable Names** | You can say it out loud |
| **Use Searchable Names** | Easy to find with Ctrl+F |
| **Avoid Encodings** | No Hungarian notation |
| **Be Consistent** | Same word for same concept |

**Remember:** Code is read far more often than it is written. Invest time in good names!
