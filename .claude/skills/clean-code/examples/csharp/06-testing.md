# C# Unit Testing Examples

> **Domain: Social Media Platform**
> Examples use User, Post, Comment, Like, and Follow entities.

## One Assert Per Test

**❌ BAD: Multiple assertions testing different behaviors**
```csharp
[Fact]
public void TestPost()
{
    var author = new User("alice", "alice@social.com");
    var post = new Post(author, "Hello World!");
    
    Assert.Equal("Hello World!", post.Content);
    Assert.Equal(author, post.Author);
    Assert.True(post.IsVisible);
    Assert.Equal(0, post.LikeCount);
    Assert.Empty(post.Comments);
}
```

**✅ GOOD: Focused tests with single responsibility**
```csharp
[Fact]
public void Post_ShouldStoreContentCorrectly()
{
    var author = new User("alice", "alice@social.com");
    var post = new Post(author, "Hello World!");
    
    Assert.Equal("Hello World!", post.Content);
}

[Fact]
public void Post_ShouldAssociateWithAuthor()
{
    var author = new User("alice", "alice@social.com");
    var post = new Post(author, "Hello World!");
    
    Assert.Equal(author.Id, post.AuthorId);
}

[Fact]
public void NewPost_ShouldBeVisibleByDefault()
{
    var post = CreatePost();
    
    Assert.True(post.IsVisible);
}

[Fact]
public void NewPost_ShouldHaveZeroLikes()
{
    var post = CreatePost();
    
    Assert.Equal(0, post.LikeCount);
}
```

## Arrange-Act-Assert Pattern

**✅ GOOD: Clear AAA structure**
```csharp
[Fact]
public void AddComment_ShouldIncreaseCommentCount()
{
    // Arrange
    var author = new User("alice", "alice@social.com");
    var commenter = new User("bob", "bob@social.com");
    var post = new Post(author, "Check out this photo!");

    // Act
    post.AddComment(commenter, "Great shot!");

    // Assert
    Assert.Equal(1, post.CommentCount);
}

[Fact]
public void LikePost_WhenNotAlreadyLiked_ShouldIncrementLikeCount()
{
    // Arrange
    var post = new PostBuilder().Build();
    var user = new UserBuilder().Build();

    // Act
    var result = post.Like(user);

    // Assert
    Assert.True(result.IsSuccess);
    Assert.Equal(1, post.LikeCount);
}

[Fact]
public void LikePost_WhenAlreadyLiked_ShouldReturnError()
{
    // Arrange
    var post = new PostBuilder().Build();
    var user = new UserBuilder().Build();
    post.Like(user); // First like

    // Act
    var result = post.Like(user); // Second like attempt

    // Assert
    Assert.True(result.IsFailure);
    Assert.Equal("User has already liked this post", result.Error);
}
```

## Use Builder Pattern for Test Data

**✅ GOOD: Fluent builders for complex test objects**
```csharp
public class UserBuilder
{
    private string _username = "testuser";
    private string _email = "test@social.com";
    private bool _isVerified = false;
    private bool _isPrivate = false;
    private int _followerCount = 0;

    public UserBuilder WithUsername(string username)
    {
        _username = username;
        return this;
    }

    public UserBuilder WithEmail(string email)
    {
        _email = email;
        return this;
    }

    public UserBuilder AsVerified()
    {
        _isVerified = true;
        return this;
    }

    public UserBuilder AsPrivateAccount()
    {
        _isPrivate = true;
        return this;
    }

    public UserBuilder WithFollowers(int count)
    {
        _followerCount = count;
        return this;
    }

    public User Build() => new(_username, _email)
    {
        IsVerified = _isVerified,
        IsPrivate = _isPrivate,
        FollowerCount = _followerCount
    };
}

public class PostBuilder
{
    private User _author = new UserBuilder().Build();
    private string _content = "Default post content";
    private List<string> _hashtags = new();
    private PostVisibility _visibility = PostVisibility.Public;
    private DateTime _createdAt = DateTime.UtcNow;

    public PostBuilder WithAuthor(User author)
    {
        _author = author;
        return this;
    }

    public PostBuilder WithContent(string content)
    {
        _content = content;
        return this;
    }

    public PostBuilder WithHashtags(params string[] hashtags)
    {
        _hashtags.AddRange(hashtags);
        return this;
    }

    public PostBuilder AsPrivate()
    {
        _visibility = PostVisibility.Private;
        return this;
    }

    public PostBuilder CreatedDaysAgo(int days)
    {
        _createdAt = DateTime.UtcNow.AddDays(-days);
        return this;
    }

    public Post Build() => new(_author, _content)
    {
        Hashtags = _hashtags,
        Visibility = _visibility,
        CreatedAt = _createdAt
    };
}
```

**Usage in tests:**
```csharp
[Fact]
public void VerifiedUser_PostsShouldShowVerificationBadge()
{
    // Arrange
    var verifiedInfluencer = new UserBuilder()
        .WithUsername("influencer")
        .AsVerified()
        .WithFollowers(100_000)
        .Build();

    var post = new PostBuilder()
        .WithAuthor(verifiedInfluencer)
        .WithContent("New product review!")
        .WithHashtags("#sponsored", "#review")
        .Build();

    // Act
    var displayInfo = post.GetDisplayInfo();

    // Assert
    Assert.True(displayInfo.ShowVerificationBadge);
}

[Fact]
public void PrivateAccount_PostsShouldOnlyBeVisibleToFollowers()
{
    // Arrange
    var privateUser = new UserBuilder()
        .AsPrivateAccount()
        .Build();

    var stranger = new UserBuilder()
        .WithUsername("stranger")
        .Build();

    var post = new PostBuilder()
        .WithAuthor(privateUser)
        .Build();

    // Act
    var canView = post.CanBeViewedBy(stranger);

    // Assert
    Assert.False(canView);
}
```

## Test Naming Convention

**✅ GOOD: MethodName_Scenario_ExpectedBehavior**
```csharp
// User actions
[Fact]
public void Follow_WhenTargetIsPublicAccount_ShouldSucceedImmediately()

[Fact]
public void Follow_WhenTargetIsPrivateAccount_ShouldCreatePendingRequest()

[Fact]
public void Follow_WhenAlreadyFollowing_ShouldReturnAlreadyFollowingError()

// Post operations
[Fact]
public void CreatePost_WithProhibitedContent_ShouldRejectWithViolationDetails()

[Fact]
public void DeletePost_ByNonAuthor_ShouldThrowUnauthorizedException()

[Fact]
public void EditPost_AfterEditWindow_ShouldReturnEditWindowExpiredError()

// Comment moderation
[Fact]
public void HideComment_ByPostAuthor_ShouldMarkAsHidden()

[Fact]
public void ReportComment_WithValidReason_ShouldQueueForModeration()
```

## Testing with Theory (Parameterized Tests)

**✅ GOOD: Using xUnit Theory for multiple scenarios**
```csharp
[Theory]
[InlineData("", false)]
[InlineData("   ", false)]
[InlineData("a", false)]                    // Too short
[InlineData("ab", false)]                   // Too short
[InlineData("abc", true)]                   // Minimum length
[InlineData("valid_username", true)]
[InlineData("user.name", true)]
[InlineData("user@name", false)]            // Invalid character
[InlineData("user name", false)]            // No spaces allowed
[InlineData("verylongusernamethatexceeds30chars", false)]
public void ValidateUsername_ShouldReturnExpectedResult(
    string username, 
    bool expectedIsValid)
{
    var result = UsernameValidator.IsValid(username);
    
    Assert.Equal(expectedIsValid, result);
}

[Theory]
[MemberData(nameof(GetHashtagTestCases))]
public void ExtractHashtags_ShouldFindAllValidHashtags(
    string content, 
    string[] expectedHashtags)
{
    var hashtags = HashtagExtractor.Extract(content);
    
    Assert.Equal(expectedHashtags, hashtags);
}

public static IEnumerable<object[]> GetHashtagTestCases()
{
    yield return new object[] { "No hashtags here", Array.Empty<string>() };
    yield return new object[] { "#hello world", new[] { "#hello" } };
    yield return new object[] { "#one #two #three", new[] { "#one", "#two", "#three" } };
    yield return new object[] { "Text #middle text", new[] { "#middle" } };
    yield return new object[] { "#CamelCase", new[] { "#CamelCase" } };
}
```

## Mocking External Dependencies

**✅ GOOD: Using Moq for notification service**
```csharp
[Fact]
public async Task CreateComment_ShouldNotifyPostAuthor()
{
    // Arrange
    var notificationServiceMock = new Mock<INotificationService>();
    var postAuthor = new UserBuilder().WithUsername("alice").Build();
    var commenter = new UserBuilder().WithUsername("bob").Build();
    var post = new PostBuilder().WithAuthor(postAuthor).Build();
    
    var commentService = new CommentService(
        notificationServiceMock.Object,
        Mock.Of<ICommentRepository>());

    // Act
    await commentService.AddComment(post.Id, commenter.Id, "Great post!");

    // Assert
    notificationServiceMock.Verify(
        x => x.SendAsync(
            postAuthor.Id,
            It.Is<Notification>(n => 
                n.Type == NotificationType.NewComment &&
                n.Message.Contains("bob"))),
        Times.Once);
}

[Fact]
public async Task MentionUser_ShouldSendMentionNotification()
{
    // Arrange
    var notificationServiceMock = new Mock<INotificationService>();
    var mentionedUser = new UserBuilder().WithUsername("charlie").Build();
    
    var userRepositoryMock = new Mock<IUserRepository>();
    userRepositoryMock
        .Setup(x => x.FindByUsernameAsync("charlie"))
        .ReturnsAsync(mentionedUser);

    var postService = new PostService(
        notificationServiceMock.Object,
        userRepositoryMock.Object);

    // Act
    await postService.CreatePost(
        authorId: Guid.NewGuid(),
        content: "Hey @charlie check this out!");

    // Assert
    notificationServiceMock.Verify(
        x => x.SendAsync(
            mentionedUser.Id,
            It.Is<Notification>(n => n.Type == NotificationType.Mention)),
        Times.Once);
}
```

## F.I.R.S.T. Principles

| Principle | Description | Example |
|-----------|-------------|---------|
| **F**ast | Tests run in milliseconds | Use in-memory repos, avoid I/O |
| **I**ndependent | No shared state between tests | Each test creates its own data |
| **R**epeatable | Same result every run | No dependency on time/random |
| **S**elf-validating | Clear pass/fail | Assert with descriptive messages |
| **T**imely | Written with production code | TDD or immediate test coverage |

**Example: Making tests repeatable (avoiding time dependency)**
```csharp
// ❌ BAD: Depends on current time
[Fact]
public void Post_ShouldBeConsideredNewIfCreatedWithin24Hours()
{
    var post = new Post(author, "content"); // Uses DateTime.Now internally
    Assert.True(post.IsNew);
}

// ✅ GOOD: Inject time abstraction
[Fact]
public void Post_CreatedWithin24Hours_ShouldBeConsideredNew()
{
    var fixedTime = new DateTime(2024, 1, 15, 12, 0, 0, DateTimeKind.Utc);
    var timeProvider = new FakeTimeProvider(fixedTime);
    
    var post = new Post(author, "content", timeProvider);
    
    Assert.True(post.IsNew); // Always returns same result
}
```

---

For integration testing strategies and architecture-level testing patterns, see:
- [Clean Architecture Testing](../../../clean-architecture/CHECKLIST.md)
- [Complete Banking Example](08-complete-example.md) for Result pattern testing
